#!/usr/bin/env python3
"""
Fetch RMS-specific expression scores from the OpenPedCan v15 gene-expression
matrix (https://github.com/d3b-center/OpenPedCan-analysis) and write a
Phase-3-compatible expression summary TSV that fills the 0.15 expression
weight in the Phase 5 confidence formula.

Source files (cached to /tmp/openpedcan/):
    gene-expression-rsem-tpm-collapsed.rds   (~285 MB, RSEM TPM matrix)
    histologies.tsv                           (~21 MB, per-sample disease labels)

Output: assets/openpedcan_expression_summary.tsv
    gene, n_rms, mean_log2tpm_rms, n_other, mean_log2tpm_other,
    log2fc_rms_vs_other, expression_zscore, expression_score, provenance

Method:
    Group OpenPedCan biospecimens with RNA-seq into RMS vs other pediatric
    cancers using the cancer_group column. For each of our 21 target genes,
    compute mean log2(TPM+1) in each group, log2 fold-change, and a
    z-score-style ratio (log2fc / pooled stdev). Map the z-score to a
    bounded expression_score in [0, 1] anchored at z=0 -> 0.0, z>=2 -> 1.0.

Dependencies:
    pyreadr (third-party; install via venv since this is a one-shot fetcher
    rather than core pipeline code). The fetcher emits a tiny TSV that the
    pipeline reads with stdlib only, so the dependency does not propagate.

Usage:
    python3 -m venv /tmp/openpedcan_venv
    /tmp/openpedcan_venv/bin/pip install pyreadr
    /tmp/openpedcan_venv/bin/python bin/fetch_openpedcan_expression.py
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
CACHE_DIR = Path("/tmp/openpedcan")
EXPR_PATH = CACHE_DIR / "gene-expression-rsem-tpm-collapsed.rds"
HIST_PATH = CACHE_DIR / "histologies.tsv"
RELEASE = "OpenPedCan v15"
RMS_CANCER_GROUP = "Rhabdomyosarcoma"

OUT_COLS = [
    "gene", "n_rms", "mean_log2tpm_rms", "n_other", "mean_log2tpm_other",
    "log2fc_rms_vs_other", "expression_zscore", "expression_score", "provenance",
]


def load_target_genes(path: Path) -> set[str]:
    out: set[str] = set()
    with path.open() as fh:
        reader = csv.DictReader((ln for ln in fh if not ln.startswith("##")), delimiter="\t")
        for row in reader:
            out.add(row["gene"].strip())
    return out


def load_histologies(path: Path) -> dict[str, str]:
    """biospecimen_id -> cancer_group label.

    Restrict to RNA-seq biospecimens since the expression matrix only covers
    those. Other rows (WGS, methylation, etc.) are dropped.
    """
    out: dict[str, str] = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("experimental_strategy", "").lower() != "rna-seq":
                continue
            bid = row.get("Kids_First_Biospecimen_ID", "").strip()
            cg = row.get("cancer_group", "").strip()
            if bid and cg:
                out[bid] = cg
    return out


def zscore_to_score(z: float) -> float:
    """Map RMS-vs-other z-score to a bounded [0, 1] expression score.

    Anchored: z <= 0 -> 0.0 (gene is NOT up-regulated in RMS), z >= 2 -> 1.0
    (gene is strongly up-regulated). Linear in between.
    """
    if z <= 0:
        return 0.0
    return min(1.0, z / 2.0)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets-kb", type=Path,
                    default=REPO_ROOT / "assets" / "targets_kb.tsv")
    ap.add_argument("--histologies", type=Path, default=HIST_PATH)
    ap.add_argument("--expression", type=Path, default=EXPR_PATH)
    ap.add_argument("--out", type=Path,
                    default=REPO_ROOT / "assets" / "openpedcan_expression_summary.tsv")
    args = ap.parse_args()

    try:
        import pyreadr  # noqa: F401
    except ImportError:
        print(
            "pyreadr not installed. This fetcher needs it to read the OpenPedCan .rds.\n"
            "Install with:\n"
            "  python3 -m venv /tmp/openpedcan_venv\n"
            "  /tmp/openpedcan_venv/bin/pip install pyreadr\n"
            "  /tmp/openpedcan_venv/bin/python bin/fetch_openpedcan_expression.py",
            file=sys.stderr,
        )
        return 2

    if not args.histologies.exists() or not args.expression.exists():
        print(
            f"missing cache files. Download with:\n"
            f"  mkdir -p {CACHE_DIR} && cd {CACHE_DIR}\n"
            f"  curl -sLO https://s3.amazonaws.com/d3b-openaccess-us-east-1-prd-pbta/open-targets/v15/histologies.tsv\n"
            f"  curl -sLO https://s3.amazonaws.com/d3b-openaccess-us-east-1-prd-pbta/open-targets/v15/gene-expression-rsem-tpm-collapsed.rds",
            file=sys.stderr,
        )
        return 1

    targets = load_target_genes(args.targets_kb)
    print(f"loaded {len(targets)} target genes", file=sys.stderr)

    bid_to_cg = load_histologies(args.histologies)
    rms_bids = {bid for bid, cg in bid_to_cg.items() if cg == RMS_CANCER_GROUP}
    other_bids = {bid for bid, cg in bid_to_cg.items() if cg and cg != RMS_CANCER_GROUP}
    print(f"RNA-seq biospecimens: {len(bid_to_cg)} total, {len(rms_bids)} RMS, {len(other_bids)} other", file=sys.stderr)

    print(f"reading {args.expression} (this takes a minute) ...", file=sys.stderr)
    import pyreadr
    rds = pyreadr.read_r(str(args.expression))
    # The .rds is a single object; pyreadr returns {None: dataframe}
    expr_df = next(iter(rds.values()))
    print(f"expression matrix shape: {expr_df.shape}", file=sys.stderr)
    print(f"first 5 columns: {list(expr_df.columns)[:5]}", file=sys.stderr)
    print(f"first 5 rows index: {list(expr_df.index)[:5]}", file=sys.stderr)

    # OpenPedCan convention: rows = gene symbols, columns = biospecimen IDs.
    # Confirm by checking if our target genes are in the index.
    expr_genes = set(expr_df.index)
    found = targets & expr_genes
    print(f"targets present in matrix: {len(found)} of {len(targets)}", file=sys.stderr)
    missing = sorted(targets - found)
    if missing:
        print(f"  missing: {missing}", file=sys.stderr)

    # Restrict to columns that have a known biospecimen mapping
    expr_cols = set(expr_df.columns)
    rms_cols = sorted(rms_bids & expr_cols)
    other_cols = sorted(other_bids & expr_cols)
    print(f"expression columns matched: {len(rms_cols)} RMS, {len(other_cols)} other", file=sys.stderr)

    rows: list[dict] = []
    for gene in sorted(targets):
        if gene not in expr_genes:
            rows.append({c: "" for c in OUT_COLS})
            rows[-1]["gene"] = gene
            rows[-1]["expression_score"] = "0.000"
            rows[-1]["provenance"] = "gene not in OpenPedCan matrix"
            continue
        rms_tpm = expr_df.loc[gene, rms_cols].astype(float).tolist() if rms_cols else []
        other_tpm = expr_df.loc[gene, other_cols].astype(float).tolist() if other_cols else []

        log_rms = [math.log2(v + 1) for v in rms_tpm]
        log_other = [math.log2(v + 1) for v in other_tpm]

        if not log_rms or not log_other:
            rows.append({
                "gene": gene, "n_rms": len(log_rms), "mean_log2tpm_rms": "",
                "n_other": len(log_other), "mean_log2tpm_other": "",
                "log2fc_rms_vs_other": "", "expression_zscore": "",
                "expression_score": "0.000",
                "provenance": "insufficient data",
            })
            continue

        m_rms = sum(log_rms) / len(log_rms)
        m_other = sum(log_other) / len(log_other)
        log2fc = m_rms - m_other
        # Pooled stdev across both groups
        all_vals = log_rms + log_other
        all_mean = sum(all_vals) / len(all_vals)
        var = sum((v - all_mean) ** 2 for v in all_vals) / max(1, len(all_vals) - 1)
        sd = math.sqrt(var) if var > 0 else 1.0
        z = log2fc / sd
        score = zscore_to_score(z)

        rows.append({
            "gene": gene,
            "n_rms": len(log_rms),
            "mean_log2tpm_rms": f"{m_rms:.3f}",
            "n_other": len(log_other),
            "mean_log2tpm_other": f"{m_other:.3f}",
            "log2fc_rms_vs_other": f"{log2fc:+.3f}",
            "expression_zscore": f"{z:+.3f}",
            "expression_score": f"{score:.3f}",
            "provenance": f"OpenPedCan_v15_RNAseq",
        })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    with args.out.open("w", newline="") as fh:
        fh.write(f"## RMS-specific expression summary (auto-generated by bin/fetch_openpedcan_expression.py)\n")
        fh.write(f"## Source: {RELEASE} gene-expression-rsem-tpm-collapsed.rds + histologies.tsv\n")
        fh.write(f"## Generated: {timestamp}\n")
        fh.write(f"## RMS RNA-seq samples: {len(rms_cols)}\n")
        fh.write(f"## Non-RMS pediatric RNA-seq samples: {len(other_cols)}\n")
        fh.write(f"## Schema: gene\\tn_rms\\tmean_log2tpm_rms\\tn_other\\tmean_log2tpm_other\\tlog2fc_rms_vs_other\\texpression_zscore\\texpression_score\\tprovenance\n")
        fh.write(f"## Method: per-gene log2(TPM+1) RMS-vs-other; z = log2fc / pooled-stdev; score = clip(z/2, 0, 1).\n")
        fh.write(f"## Higher expression_score => gene is more up-regulated in RMS vs other pediatric cancers.\n")
        fh.write(f"## To regenerate: see fetcher docstring (requires pyreadr in a venv).\n")
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
