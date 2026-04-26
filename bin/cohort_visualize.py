#!/usr/bin/env python3
"""
Cohort-level visualizations for the RMS-ISP TARGET-RT runner.

Reads results/target_rt_cohort_summary.tsv plus per-sample p5.tsv outputs,
writes a long-format aggregation TSV and three SVG charts:

  - cohort_mechanisms.svg     horizontal bar chart of top_mechanism counts
  - cohort_druggability.svg   gene x subtype fraction matrix (always 21 x 4)
  - cohort_per_sample.svg     gene x sample heatmap (only when N <= 100)

Pure stdlib. Called from bin/run_target_rt.py at the end of a cohort run, or
runnable standalone for fast viz iteration without re-running 36 pipelines.
"""
from __future__ import annotations

import csv
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# Confidence threshold: cells below this count as "not druggable" and render blank.
# Pinned to 0.10 to match the passenger-sanity assertion in tests/cases.toml.
DRUGGABILITY_THRESHOLD = 0.10

# Color choices (committed here so the test can assert on hex values).
HEATMAP_RAMP = ["#eff3ff", "#bdd7e7", "#6baed6", "#3182bd", "#08519c"]
HEATMAP_BINS = [(0.10, 0.20), (0.20, 0.40), (0.40, 0.60), (0.60, 0.80), (0.80, 1.01)]
EMPTY_CELL = "#f5f5f5"
SUBTYPE_PALETTE = {"FN": "#377eb8", "FP": "#e41a1c", "ALL": "#984ea3"}

# Suppress per-sample heatmap above this cohort size (illegible).
PER_SAMPLE_MAX_N = 100


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def _read_cohort_summary(path: Path) -> list[dict]:
    """Returns the per-sample metadata rows from the cohort summary TSV."""
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def _read_p5(path: Path) -> list[dict]:
    """Returns all event x drug rows from a phase-5 TSV."""
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def aggregate_gene_matrix(cohort_tsv: Path, target_rt_dir: Path) -> list[dict]:
    """
    Returns long-format rows: {sample_id, study, subtype, gene, max_confidence}.

    Sample list comes from the cohort summary TSV (authoritative; samples that
    produced no scored rows are correctly excluded). For each sample, reads the
    paired p5.tsv and groups its rows by (sample, gene), keeping the max
    confidence. Cells below DRUGGABILITY_THRESHOLD are dropped.

    Output is sorted (sample_id, gene) for determinism.
    """
    cohort = _read_cohort_summary(cohort_tsv)
    out: list[dict] = []
    for sample_meta in cohort:
        sid = sample_meta["sample_id"]
        p5_path = target_rt_dir / sid / "p5.tsv"
        if not p5_path.exists():
            continue
        per_gene: dict[str, float] = {}
        for r in _read_p5(p5_path):
            gene = r.get("gene", "")
            try:
                conf = float(r.get("confidence", "") or 0)
            except ValueError:
                conf = 0.0
            if not gene or conf < DRUGGABILITY_THRESHOLD:
                continue
            if conf > per_gene.get(gene, 0.0):
                per_gene[gene] = conf
        for gene, conf in per_gene.items():
            out.append({
                "sample_id": sid,
                "study": sample_meta.get("study", ""),
                "subtype": sample_meta.get("subtype", ""),
                "gene": gene,
                "max_confidence": f"{conf:.3f}",
            })
    out.sort(key=lambda r: (r["sample_id"], r["gene"]))
    return out


def write_gene_matrix_tsv(rows: list[dict], path: Path) -> None:
    """Writes the aggregation rows to a tab-separated file with a fixed column order."""
    cols = ["sample_id", "study", "subtype", "gene", "max_confidence"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})


# ---------------------------------------------------------------------------
# SVG primitives
# ---------------------------------------------------------------------------

def _svg_root(width: int, height: int) -> ET.Element:
    """Returns an <svg> element with sensible defaults."""
    return ET.Element("svg", {
        "xmlns": "http://www.w3.org/2000/svg",
        "viewBox": f"0 0 {width} {height}",
        "width": str(width),
        "height": str(height),
        "font-family": "system-ui, -apple-system, sans-serif",
        "font-size": "12",
    })


# ---------------------------------------------------------------------------
# Mechanism distribution
# ---------------------------------------------------------------------------

def render_mechanism_chart(mechanism_counts: dict[str, int], cohort_size: int) -> str:
    """
    Horizontal bar chart of top mechanism counts.

    Mechanisms with count < max(1, cohort_size // 100) are binned into a single
    'other (N items)' bar at the bottom (omitted if no such items).
    """
    items = sorted(mechanism_counts.items(), key=lambda kv: (-kv[1], kv[0]))
    threshold = max(1, cohort_size // 100)
    head = [(m, c) for m, c in items if c >= threshold]
    tail = [(m, c) for m, c in items if c < threshold]
    bars = list(head)
    if tail:
        bars.append((f"other ({len(tail)} items)", sum(c for _, c in tail)))

    bar_h = 22
    pad_top = 30
    pad_bottom = 20
    pad_left = 220
    pad_right = 50
    plot_w = 400
    width = pad_left + plot_w + pad_right
    height = pad_top + pad_bottom + max(1, len(bars)) * bar_h

    svg = _svg_root(width, height)
    title = ET.SubElement(svg, "text", {
        "x": str(width // 2), "y": "20",
        "text-anchor": "middle", "font-weight": "bold",
    })
    title.text = f"Top mechanism per sample (N={cohort_size})"

    if not bars:
        empty = ET.SubElement(svg, "text", {
            "x": str(width // 2), "y": str(height // 2),
            "text-anchor": "middle",
        })
        empty.text = "(no mechanisms)"
        return ET.tostring(svg, encoding="unicode")

    max_count = max(c for _, c in bars)
    for i, (mech, count) in enumerate(bars):
        y = pad_top + i * bar_h
        bar_w = int(plot_w * (count / max_count))
        ET.SubElement(svg, "text", {
            "x": str(pad_left - 6), "y": str(y + bar_h - 7),
            "text-anchor": "end",
        }).text = mech
        ET.SubElement(svg, "rect", {
            "class": "bar",
            "data-label": mech,
            "x": str(pad_left), "y": str(y + 3),
            "width": str(bar_w), "height": str(bar_h - 6),
            "fill": "#3182bd",
        })
        ET.SubElement(svg, "text", {
            "x": str(pad_left + bar_w + 4), "y": str(y + bar_h - 7),
        }).text = str(count)

    return ET.tostring(svg, encoding="unicode")


# ---------------------------------------------------------------------------
# Per-target druggability (gene x subtype, always 21 x 4)
# ---------------------------------------------------------------------------

def _load_target_genes(path: Path) -> list[str]:
    """Reads gene names from assets/targets_kb.tsv, skipping ## header lines."""
    out: list[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith("##") or line.strip() == "":
                continue
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            head = parts[0]
            if head == "gene":
                continue
            out.append(head)
    return out


def compute_druggability_matrix(
    agg_rows: list[dict],
    cohort_meta: list[dict],
    genes: list[str],
) -> dict[tuple[str, str], float]:
    """
    For each (gene, subtype-bucket), returns the fraction of samples in that
    bucket with a druggable hit on that gene.

    Buckets: 'FN', 'FP', 'ALL' (samples whose recorded subtype is exactly
    'ALL'), and 'whole_cohort' (every sample in cohort_meta).
    """
    buckets = ["FN", "FP", "ALL", "whole_cohort"]

    sample_buckets: dict[str, set[str]] = {b: set() for b in buckets}
    for s in cohort_meta:
        sid = s["sample_id"]
        sub = s.get("subtype", "")
        if sub in ("FN", "FP", "ALL"):
            sample_buckets[sub].add(sid)
        sample_buckets["whole_cohort"].add(sid)

    by_gene: dict[str, set[str]] = {}
    for r in agg_rows:
        by_gene.setdefault(r["gene"], set()).add(r["sample_id"])

    matrix: dict[tuple[str, str], float] = {}
    for gene in genes:
        hit_samples = by_gene.get(gene, set())
        for b in buckets:
            denom = len(sample_buckets[b])
            if denom == 0:
                matrix[(gene, b)] = 0.0
                continue
            num = len(hit_samples & sample_buckets[b])
            matrix[(gene, b)] = num / denom
    return matrix


def _color_for_value(value: float) -> str:
    """Maps a 0..1 fraction to the heatmap ramp; below threshold -> EMPTY_CELL."""
    if value < DRUGGABILITY_THRESHOLD:
        return EMPTY_CELL
    for (lo, hi), color in zip(HEATMAP_BINS, HEATMAP_RAMP):
        if lo <= value < hi:
            return color
    return HEATMAP_RAMP[-1]


def render_druggability_chart(
    matrix: dict[tuple[str, str], float],
    genes: list[str],
    subtype_columns: list[str],
) -> str:
    """
    Gene rows x subtype columns. Always 21 x 4 in production. Cell color is the
    fraction; cell text is the fraction to two decimals when above threshold.
    """
    cell_w = 90
    cell_h = 22
    pad_top = 60
    pad_bottom = 20
    pad_left = 90
    pad_right = 20
    width = pad_left + cell_w * len(subtype_columns) + pad_right
    height = pad_top + cell_h * len(genes) + pad_bottom

    svg = _svg_root(width, height)
    ET.SubElement(svg, "text", {
        "x": str(width // 2), "y": "22",
        "text-anchor": "middle", "font-weight": "bold",
    }).text = "Per-target cohort druggability (fraction of subtype with confidence >= 0.10)"

    for j, sub in enumerate(subtype_columns):
        cx = pad_left + j * cell_w + cell_w // 2
        ET.SubElement(svg, "text", {
            "x": str(cx), "y": str(pad_top - 8),
            "text-anchor": "middle", "font-weight": "bold",
        }).text = sub

    for i, gene in enumerate(genes):
        ry = pad_top + i * cell_h
        ET.SubElement(svg, "text", {
            "x": str(pad_left - 6), "y": str(ry + cell_h - 7),
            "text-anchor": "end",
        }).text = gene
        for j, sub in enumerate(subtype_columns):
            cx = pad_left + j * cell_w
            value = matrix.get((gene, sub), 0.0)
            ET.SubElement(svg, "rect", {
                "class": "cell",
                "data-gene": gene,
                "data-subtype": sub,
                "x": str(cx), "y": str(ry),
                "width": str(cell_w), "height": str(cell_h),
                "fill": _color_for_value(value),
                "stroke": "#ffffff", "stroke-width": "1",
            })
            if value >= DRUGGABILITY_THRESHOLD:
                ET.SubElement(svg, "text", {
                    "x": str(cx + cell_w // 2),
                    "y": str(ry + cell_h - 7),
                    "text-anchor": "middle",
                    "fill": "#ffffff" if value >= 0.40 else "#000000",
                }).text = f"{value:.2f}"
    return ET.tostring(svg, encoding="unicode")


# ---------------------------------------------------------------------------
# Per-sample heatmap (conditional on N <= PER_SAMPLE_MAX_N)
# ---------------------------------------------------------------------------

def render_per_sample_heatmap(
    agg_rows: list[dict],
    cohort_meta: list[dict],
    genes: list[str],
) -> str | None:
    """
    Gene rows x sample columns, cell color = max_confidence.

    Returns None if len(cohort_meta) > PER_SAMPLE_MAX_N (the chart would be
    illegible). Caller is expected to handle None by writing a stub note in
    the surrounding markdown instead of an <img> tag.

    Sample columns are sorted by (subtype, sample_id). Sample-ID labels
    appear only when len(cohort_meta) <= 50.
    """
    if len(cohort_meta) > PER_SAMPLE_MAX_N:
        return None

    samples = sorted(cohort_meta, key=lambda s: (s.get("subtype", ""), s["sample_id"]))
    sids = [s["sample_id"] for s in samples]
    show_labels = len(samples) <= 50

    by_key: dict[tuple[str, str], float] = {}
    for r in agg_rows:
        try:
            by_key[(r["sample_id"], r["gene"])] = float(r["max_confidence"])
        except (KeyError, ValueError):
            continue

    cell_w = 16
    cell_h = 18
    pad_top = 70
    pad_bottom = 60 if show_labels else 30
    pad_left = 90
    pad_right = 10
    stripe_h = 10
    width = pad_left + cell_w * len(samples) + pad_right
    height = pad_top + cell_h * len(genes) + pad_bottom

    svg = _svg_root(width, height)
    ET.SubElement(svg, "text", {
        "x": str(width // 2), "y": "22",
        "text-anchor": "middle", "font-weight": "bold",
    }).text = f"Per-sample druggability (N={len(samples)})"

    stripe_y = pad_top - stripe_h - 4
    for j, s in enumerate(samples):
        x = pad_left + j * cell_w
        ET.SubElement(svg, "rect", {
            "class": "subtype-stripe",
            "data-subtype": s.get("subtype", ""),
            "x": str(x), "y": str(stripe_y),
            "width": str(cell_w), "height": str(stripe_h),
            "fill": SUBTYPE_PALETTE.get(s.get("subtype", ""), "#cccccc"),
        })

    for i, gene in enumerate(genes):
        ry = pad_top + i * cell_h
        ET.SubElement(svg, "text", {
            "x": str(pad_left - 6), "y": str(ry + cell_h - 5),
            "text-anchor": "end",
        }).text = gene
        for j, sid in enumerate(sids):
            cx = pad_left + j * cell_w
            value = by_key.get((sid, gene), 0.0)
            ET.SubElement(svg, "rect", {
                "class": "cell",
                "data-sample": sid,
                "data-gene": gene,
                "x": str(cx), "y": str(ry),
                "width": str(cell_w), "height": str(cell_h),
                "fill": _color_for_value(value),
                "stroke": "#ffffff", "stroke-width": "0.5",
            })

    if show_labels:
        label_y = pad_top + cell_h * len(genes) + 12
        for j, sid in enumerate(sids):
            cx = pad_left + j * cell_w + cell_w // 2
            t = ET.SubElement(svg, "text", {
                "x": str(cx), "y": str(label_y),
                "text-anchor": "end",
                "transform": f"rotate(-60 {cx} {label_y})",
                "font-size": "10",
            })
            t.text = sid

    return ET.tostring(svg, encoding="unicode")


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def main(
    cohort_tsv: Path,
    target_rt_dir: Path,
    out_dir: Path,
    targets_kb: Path,
) -> dict:
    """
    Generate the aggregation TSV and the three SVG files.

    Returns a status dict consumed by bin/run_target_rt.py:
        {
          "n_samples": int,
          "mechanisms": Path,
          "druggability": Path,
          "per_sample": Path | None,
        }
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    cohort_meta = _read_cohort_summary(cohort_tsv)

    agg_rows = aggregate_gene_matrix(cohort_tsv, target_rt_dir)
    matrix_path = out_dir / "cohort_gene_matrix.tsv"
    write_gene_matrix_tsv(agg_rows, matrix_path)

    mech_counts = Counter(s["top_mechanism"] for s in cohort_meta
                          if s.get("top_mechanism"))
    mech_path = out_dir / "cohort_mechanisms.svg"
    mech_path.write_text(render_mechanism_chart(dict(mech_counts),
                                                cohort_size=len(cohort_meta)))

    genes = _load_target_genes(targets_kb)
    drug_matrix = compute_druggability_matrix(agg_rows, cohort_meta, genes)
    drug_path = out_dir / "cohort_druggability.svg"
    drug_path.write_text(render_druggability_chart(
        drug_matrix, genes, ["FN", "FP", "ALL", "whole_cohort"]))

    per_sample_svg = render_per_sample_heatmap(agg_rows, cohort_meta, genes)
    per_sample_path: Path | None = None
    if per_sample_svg is not None:
        per_sample_path = out_dir / "cohort_per_sample.svg"
        per_sample_path.write_text(per_sample_svg)

    return {
        "n_samples": len(cohort_meta),
        "mechanisms": mech_path,
        "druggability": drug_path,
        "per_sample": per_sample_path,
    }


def _cli() -> int:
    import argparse
    import sys
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cohort-tsv", type=Path,
                    default=REPO_ROOT / "results" / "target_rt_cohort_summary.tsv")
    ap.add_argument("--target-rt-dir", type=Path,
                    default=REPO_ROOT / "results" / "target_rt")
    ap.add_argument("--out-dir", type=Path,
                    default=REPO_ROOT / "results" / "target_rt")
    ap.add_argument("--targets-kb", type=Path,
                    default=REPO_ROOT / "assets" / "targets_kb.tsv")
    args = ap.parse_args()
    if not args.cohort_tsv.exists():
        print(f"cohort summary not found at {args.cohort_tsv}; "
              f"run bin/run_target_rt.py first", file=sys.stderr)
        return 1
    status = main(args.cohort_tsv, args.target_rt_dir, args.out_dir,
                  args.targets_kb)
    print(f"wrote: {status['mechanisms']}")
    print(f"wrote: {status['druggability']}")
    if status["per_sample"]:
        print(f"wrote: {status['per_sample']}")
    else:
        print(f"per-sample heatmap suppressed at N={status['n_samples']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_cli())
