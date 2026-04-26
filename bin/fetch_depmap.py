#!/usr/bin/env python3
"""
Fetch real DepMap CRISPR Chronos gene-effect data for the RMS cell-line panel
and the RMS-ISP target gene set, and write a Phase-3-compatible summary TSV.

Source: DepMap Public 26Q1 Chronos parameters published on Figshare
(article 31660582). The gene_effect.csv file (~431 MB) holds per-cell-line
per-gene Chronos scores. We stream-download it (with caching) and project to
just the (RMS rows × target columns) subset we care about.

Cell-line registry: assets/rms_cell_lines.tsv (curated, with FP/FN subtype).
Target gene registry: assets/targets_kb.tsv (the 21 RMS-relevant genes).

Output schema matches assets/depmap_rms_summary.tsv exactly so Phase 3 can
read the result without any code change. The script REPLACES that file in
place when run with --write-summary.

Usage:
    bin/fetch_depmap.py                          # download + summarize, write to default location
    bin/fetch_depmap.py --no-download            # use cached /tmp/depmap/gene_effect.csv if present
    bin/fetch_depmap.py --out /tmp/peek.tsv      # write somewhere else for inspection

Network: one HTTP GET to figshare.com (~431 MB). Cached to /tmp/depmap/ so
re-runs are free. The bundled summary TSV is small and committed to git.
"""
from __future__ import annotations

import argparse
import csv
import statistics
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
DEPMAP_URL = "https://ndownloader.figshare.com/files/62677015"
DEPMAP_RELEASE = "DepMap Public 26Q1 (Figshare article 31660582)"
CACHE_PATH = Path("/tmp/depmap/gene_effect.csv")
ESSENTIAL_THRESHOLD = -0.5  # Chronos < -0.5 is the standard "essential" cut-off

OUT_COLS = [
    "gene", "n_lines", "mean_chronos_all", "mean_chronos_fp",
    "mean_chronos_fn", "pct_essential", "dependency_score", "provenance",
]


def download_if_needed(url: str, dest: Path, force: bool) -> Path:
    if dest.exists() and not force:
        size_mb = dest.stat().st_size / 1e6
        print(f"using cached {dest} ({size_mb:.1f} MB)", file=sys.stderr)
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"downloading {url} -> {dest} ...", file=sys.stderr)
    with urllib.request.urlopen(url, timeout=600) as resp, dest.open("wb") as fh:
        total = 0
        while chunk := resp.read(1 << 20):
            fh.write(chunk)
            total += len(chunk)
        print(f"  downloaded {total / 1e6:.1f} MB", file=sys.stderr)
    return dest


def load_cell_lines(path: Path) -> dict[str, dict]:
    """ACH-ID -> {cell_line, subtype, fusion}."""
    out: dict[str, dict] = {}
    with path.open() as fh:
        reader = csv.DictReader((ln for ln in fh if not ln.startswith("##")), delimiter="\t")
        for row in reader:
            out[row["model_id"].strip()] = {
                "cell_line": row["cell_line"].strip(),
                "subtype": row["subtype"].strip(),
                "fusion": row["fusion"].strip(),
            }
    return out


def load_target_genes(path: Path) -> set[str]:
    out: set[str] = set()
    with path.open() as fh:
        reader = csv.DictReader((ln for ln in fh if not ln.startswith("##")), delimiter="\t")
        for row in reader:
            out.add(row["gene"].strip())
    return out


def parse_gene_effect(
    gene_effect_path: Path,
    cell_lines: dict[str, dict],
    target_genes: set[str],
) -> dict[str, list[tuple[float, str]]]:
    """Stream gene_effect.csv. Return {gene: [(chronos, subtype), ...]}.

    Header row: empty,GENE_SYMBOL (entrez_id),GENE_SYMBOL (entrez_id),...
    Data rows:  ACH-XXXXXX,float,float,...
    """
    chronos_by_gene: dict[str, list[tuple[float, str]]] = {g: [] for g in target_genes}
    used_models: set[str] = set()

    with gene_effect_path.open() as fh:
        header = fh.readline().rstrip("\n").split(",")
        # Map column index -> gene symbol (strip the " (entrez_id)" suffix)
        col_idx_to_gene: dict[int, str] = {}
        for i, h in enumerate(header):
            if i == 0:
                continue  # row identifier column
            sym = h.split(" (")[0].strip()
            if sym in target_genes:
                col_idx_to_gene[i] = sym

        missing_targets = target_genes - set(col_idx_to_gene.values())
        if missing_targets:
            print(f"WARN: targets absent from gene_effect.csv header: {sorted(missing_targets)}", file=sys.stderr)

        for line in fh:
            parts = line.rstrip("\n").split(",")
            ach = parts[0]
            if ach not in cell_lines:
                continue
            used_models.add(ach)
            subtype = cell_lines[ach]["subtype"]
            for col_idx, gene in col_idx_to_gene.items():
                if col_idx >= len(parts):
                    continue
                cell = parts[col_idx]
                if not cell or cell.lower() == "na":
                    continue
                try:
                    chronos_by_gene[gene].append((float(cell), subtype))
                except ValueError:
                    continue

    print(f"used {len(used_models)} of {len(cell_lines)} registered RMS lines", file=sys.stderr)
    missing = sorted(set(cell_lines) - used_models)
    if missing:
        print(f"WARN: registered RMS lines absent from gene_effect.csv: {missing}", file=sys.stderr)
    return chronos_by_gene


def summarize(
    chronos_by_gene: dict[str, list[tuple[float, str]]],
    cell_lines: dict[str, dict],
) -> list[dict]:
    """Aggregate per-gene Chronos values to the depmap_rms_summary schema."""
    n_lines_total = len(chronos_by_gene_lines := sorted(set(
        ach for ach, meta in cell_lines.items()
    )))
    rows: list[dict] = []
    for gene in sorted(chronos_by_gene):
        values = chronos_by_gene[gene]
        if not values:
            # Gene not in DepMap header at all; emit a blank row so the schema stays consistent.
            rows.append({
                "gene": gene, "n_lines": 0,
                "mean_chronos_all": "", "mean_chronos_fp": "", "mean_chronos_fn": "",
                "pct_essential": "", "dependency_score": "0.000",
                "provenance": "gene not in DepMap 26Q1",
            })
            continue
        all_vals = [v for v, _ in values]
        fp_vals = [v for v, s in values if s == "FP"]
        fn_vals = [v for v, s in values if s == "FN"]
        n = len(all_vals)
        mean_all = statistics.mean(all_vals)
        mean_fp = statistics.mean(fp_vals) if fp_vals else mean_all
        mean_fn = statistics.mean(fn_vals) if fn_vals else mean_all
        pct_essential = sum(1 for v in all_vals if v < ESSENTIAL_THRESHOLD) / n
        dependency_score = chronos_to_score(mean_all)
        rows.append({
            "gene": gene,
            "n_lines": n,
            "mean_chronos_all": f"{mean_all:.3f}",
            "mean_chronos_fp": f"{mean_fp:.3f}",
            "mean_chronos_fn": f"{mean_fn:.3f}",
            "pct_essential": f"{pct_essential:.3f}",
            "dependency_score": f"{dependency_score:.3f}",
            "provenance": "DepMap_26Q1_Chronos",
        })
    return rows


def chronos_to_score(mean_chronos: float) -> float:
    """Map Chronos (~ +0.5 to -2.0) to a 0..1 dependency score.

    Convention: Chronos > 0  -> not essential -> score near 0.
                Chronos = -0.5 -> just past essentiality threshold -> ~0.5.
                Chronos = -1.5 -> strongly essential -> ~1.0.
    Linear ramp clipped to [0, 1] anchored at those two reference points.
    """
    # Map Chronos = 0 to score 0; Chronos = -1.5 to score 1.0.
    score = max(0.0, min(1.0, -mean_chronos / 1.5))
    return score


def write_summary(rows: list[dict], out_path: Path, *, n_lines: int) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    with out_path.open("w", newline="") as fh:
        fh.write("## RMS-ISP DepMap Dependency Summary (REAL DATA, generated by bin/fetch_depmap.py)\n")
        fh.write(f"## Source: {DEPMAP_RELEASE}\n")
        fh.write(f"## File: {DEPMAP_URL}\n")
        fh.write(f"## Generated: {timestamp}\n")
        fh.write(f"## Cell-line registry: assets/rms_cell_lines.tsv ({n_lines} RMS lines used)\n")
        fh.write(f"## Essentiality threshold: Chronos < {ESSENTIAL_THRESHOLD}\n")
        fh.write(f"## Schema: gene\\tn_lines\\tmean_chronos_all\\tmean_chronos_fp\\tmean_chronos_fn\\tpct_essential\\tdependency_score\\tprovenance\n")
        fh.write("##\n")
        fh.write("## - mean_chronos_*: mean Chronos score across RMS cell lines (lower = more essential; <-0.5 = essential, <-1.0 = strongly essential)\n")
        fh.write("## - pct_essential: fraction of RMS lines where Chronos < -0.5\n")
        fh.write("## - dependency_score: linear map of mean_chronos_all to [0,1], anchored at Chronos=0 -> 0.0 and Chronos=-1.5 -> 1.0; clipped\n")
        fh.write("##\n")
        fh.write("## To regenerate: bin/fetch_depmap.py (downloads ~431 MB from Figshare; caches to /tmp/depmap/)\n")
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell-lines", type=Path,
                    default=REPO_ROOT / "assets" / "rms_cell_lines.tsv")
    ap.add_argument("--targets-kb", type=Path,
                    default=REPO_ROOT / "assets" / "targets_kb.tsv")
    ap.add_argument("--gene-effect", type=Path, default=CACHE_PATH,
                    help="Path to (cached) gene_effect.csv. Downloaded if missing.")
    ap.add_argument("--out", type=Path,
                    default=REPO_ROOT / "assets" / "depmap_rms_summary.tsv")
    ap.add_argument("--no-download", action="store_true",
                    help="Fail rather than download if cache is missing.")
    ap.add_argument("--force-download", action="store_true",
                    help="Re-download even if cached.")
    args = ap.parse_args()

    if args.no_download and not args.gene_effect.exists():
        print(f"--no-download set but {args.gene_effect} missing", file=sys.stderr)
        return 1

    if not args.no_download:
        download_if_needed(DEPMAP_URL, args.gene_effect, args.force_download)

    cell_lines = load_cell_lines(args.cell_lines)
    targets = load_target_genes(args.targets_kb)
    print(f"registry: {len(cell_lines)} RMS lines, {len(targets)} target genes", file=sys.stderr)

    chronos_by_gene = parse_gene_effect(args.gene_effect, cell_lines, targets)
    rows = summarize(chronos_by_gene, cell_lines)

    n_used = max((int(r["n_lines"]) for r in rows if r["n_lines"]), default=0)
    write_summary(rows, args.out, n_lines=n_used)
    print(f"wrote {args.out} ({len(rows)} gene rows from {n_used} RMS lines)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
