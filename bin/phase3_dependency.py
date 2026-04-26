#!/usr/bin/env python3
"""
Phase 3: Expression and dependency integration.

Joins each variant's gene against assets/depmap_rms_summary.tsv and adds the
DepMap dependency score plus the FP- vs FN-stratified Chronos values.

The bundled summary is a curated PLACEHOLDER (see assets/README.md). A real
DepMap 24Q2+ pull replaces this file with no code change.

A `--subtype` flag (FP|FN|UNKNOWN) lets the report use the subtype-appropriate
Chronos column. v0.1 defaults to ALL since the toy patient is unsubtyped.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_depmap(path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    with path.open() as fh:
        reader = csv.DictReader((ln for ln in fh if not ln.startswith("##")), delimiter="\t")
        for row in reader:
            out[row["gene"].strip()] = row
    return out


def load_expression(path: Path | None) -> dict[str, dict]:
    if not path or not path.exists() or path.stat().st_size == 0:
        return {}
    out: dict[str, dict] = {}
    with path.open() as fh:
        reader = csv.DictReader((ln for ln in fh if not ln.startswith("##")), delimiter="\t")
        for row in reader:
            out[row["gene"].strip()] = row
    return out


def _load_opentargets_lookup(efo_id: str, cache_dir: Path | None):
    """Returns a callable gene_symbol -> dict | None, with in-memory memoization."""
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).resolve().parent))
    import opentargets_client as ot
    cdir = cache_dir or ot.DEFAULT_CACHE_DIR
    memo: dict[str, dict | None] = {}

    def _lookup(symbol: str) -> dict | None:
        if symbol in memo:
            return memo[symbol]
        result = ot.lookup_gene_disease(symbol, efo_id=efo_id, cache_dir=cdir)
        memo[symbol] = result
        return result

    return _lookup


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, type=Path)
    ap.add_argument("--depmap", required=True, type=Path)
    ap.add_argument("--expression", type=Path, default=None,
                    help="Optional OpenPedCan expression summary (output of bin/fetch_openpedcan_expression.py); fills the 0.15 expression weight in Phase 5.")
    ap.add_argument("--subtype", default="ALL", choices=["ALL", "FP", "FN", "UNKNOWN"])
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--opentargets-efo-id", default="",
                    help="EFO ID for the disease (e.g., EFO_0002918 for "
                         "rhabdomyosarcoma). When set, phase 3 looks up "
                         "OpenTargets gene-disease association score and "
                         "known-drug count per gene; when empty, the columns "
                         "stay empty. v0.17 informational; not in confidence "
                         "formula yet.")
    ap.add_argument("--opentargets-cache-dir", type=Path, default=None,
                    help="Cache dir for OpenTargets responses. "
                         "Default: <repo>/data/opentargets_cache.")
    args = ap.parse_args()

    depmap = load_depmap(args.depmap)
    expression = load_expression(args.expression)

    with args.inp.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        in_cols = list(reader.fieldnames or [])

    add_cols = [
        "depmap_n_lines", "depmap_chronos_all", "depmap_chronos_fp",
        "depmap_chronos_fn", "depmap_pct_essential", "dependency_score",
        "dependency_reason",
        "expression_log2fc", "expression_zscore", "expression_score", "expression_reason",
        "opentargets_score", "opentargets_disease",
    ]
    out_cols = in_cols + add_cols

    # OpenTargets lookup is optional; lazy-imported so the phase script does
    # not pull urllib unnecessarily when the flag is empty.
    ot_lookup = _load_opentargets_lookup(
        args.opentargets_efo_id, args.opentargets_cache_dir,
    ) if args.opentargets_efo_id else None

    chronos_field = {
        "ALL": "mean_chronos_all",
        "FP": "mean_chronos_fp",
        "FN": "mean_chronos_fn",
        "UNKNOWN": "mean_chronos_all",
    }[args.subtype]

    # Oncogene-addiction floor. The DepMap RMS cell-line panel does not always
    # include lines carrying the relevant activating mutation for every oncogene
    # (e.g., the panel may have no FGFR4 V550L lines). Panel-level Chronos
    # therefore understates dependency in a tumor that DOES carry the activating
    # event. Per the oncogene-addiction literature, a confirmed activating
    # DRIVER in an oncogene creates a meaningful dependency on that gene
    # regardless of what wild-type panel members show. We encode this as a
    # minimum dependency_score for (call=DRIVER, role=oncogene) rows.
    ONCOGENE_DRIVER_FLOOR = 0.50

    for r in rows:
        gene = r.get("gene", "")
        d = depmap.get(gene)
        if d is None:
            r["depmap_n_lines"] = ""
            r["depmap_chronos_all"] = ""
            r["depmap_chronos_fp"] = ""
            r["depmap_chronos_fn"] = ""
            r["depmap_pct_essential"] = ""
            r["dependency_score"] = "0.000"
            r["dependency_reason"] = "gene not in DepMap RMS summary"
            continue
        r["depmap_n_lines"] = d["n_lines"]
        r["depmap_chronos_all"] = d["mean_chronos_all"]
        r["depmap_chronos_fp"] = d["mean_chronos_fp"]
        r["depmap_chronos_fn"] = d["mean_chronos_fn"]
        r["depmap_pct_essential"] = d["pct_essential"]
        # The dependency_score is precomputed in the summary; honour the bundled value.
        # Subtype-aware override: if the subtype-specific Chronos is meaningfully more
        # essential than ALL, scale the score upward (capped at 1.0).
        base = float(d["dependency_score"])
        chronos_used = float(d.get(chronos_field, d["mean_chronos_all"]))
        chronos_all = float(d["mean_chronos_all"])
        bonus = max(0.0, (chronos_all - chronos_used) * 0.3)  # 0.3 per Chronos unit improvement
        score = min(1.0, base + bonus)
        reason = (
            f"{d['provenance']}; Chronos[{args.subtype}]={chronos_used:+.2f}, "
            f"%essential={float(d['pct_essential']) * 100:.0f}%"
        )
        if r.get("call") == "DRIVER" and r.get("role") == "oncogene" and score < ONCOGENE_DRIVER_FLOOR:
            reason += f"; oncogene-addiction floor applied (panel score {score:.2f} -> {ONCOGENE_DRIVER_FLOOR:.2f})"
            score = ONCOGENE_DRIVER_FLOOR
        r["dependency_score"] = f"{score:.3f}"
        r["dependency_reason"] = reason

        # Expression layer (optional; only filled when --expression supplied)
        e = expression.get(gene)
        if e:
            r["expression_log2fc"] = e.get("log2fc_rms_vs_other", "")
            r["expression_zscore"] = e.get("expression_zscore", "")
            r["expression_score"] = e.get("expression_score", "0.000")
            r["expression_reason"] = (
                f"{e.get('provenance','')}; mean log2(TPM+1) RMS={e.get('mean_log2tpm_rms','')} "
                f"vs other={e.get('mean_log2tpm_other','')} (n_rms={e.get('n_rms','')})"
            )
        else:
            r["expression_log2fc"] = ""
            r["expression_zscore"] = ""
            r["expression_score"] = "0.000"
            r["expression_reason"] = "no expression data" if not expression else "gene not in expression summary"

        # OpenTargets gene-disease association (informational; not in formula).
        if ot_lookup is not None and gene:
            ot = ot_lookup(gene)
            if ot:
                r["opentargets_score"] = f"{ot['association_score']:.3f}"
                r["opentargets_disease"] = ot.get("matched_disease_name", "")
            else:
                r["opentargets_score"] = ""
                r["opentargets_disease"] = ""
        else:
            r["opentargets_score"] = ""
            r["opentargets_disease"] = ""

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
