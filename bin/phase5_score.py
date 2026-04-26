#!/usr/bin/env python3
"""
Phase 5: Confidence scoring and patient-level report generation.

Combines the Phase 1-4 component scores into a per-variant-per-drug confidence
score using the weights defined in the pilot plan §4 Aim 2 Phase 5:

    variant pathogenicity   0.25
    structural match        0.15
    dependency              0.25
    expression specificity  0.15   (deferred in v0.1; contributes 0)
    drug-matching evidence  0.20

The report explicitly exposes every component so a clinician can audit the
call. Top-ranked target-drug pairs are listed first; per-variant detail follows.

Outputs:
    {sample}.phase5.tsv   long-format scored table (one row per variant-drug)
    {sample}.report.md    human-readable markdown report
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
from datetime import datetime, timezone
from pathlib import Path


WEIGHTS = {
    "variant": 0.25,
    "structural": 0.15,
    "dependency": 0.25,
    "expression": 0.15,
    "drug": 0.20,
}


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def safe_float(s: str, default: float = 0.0) -> float:
    try:
        return float(s)
    except (TypeError, ValueError):
        return default


def event_label(r: dict) -> str:
    """Render a short human-readable label for a row's molecular event."""
    et = r.get("event_type", "")
    if et == "snv":
        return r.get("hgvsp_short") or r.get("consequence", "")
    if et == "cna":
        cons = r.get("consequence", "")
        cn = r.get("copy_number", "")
        return f"{cons} (CN={cn})" if cn else cons
    if et == "fusion":
        partner = r.get("fusion_partner", "")
        return f"fusion::{partner}" if partner else "fusion"
    return r.get("hgvsp_short") or r.get("consequence", "")


def confidence(row: dict) -> tuple[float, dict[str, float]]:
    components = {
        "variant": safe_float(row.get("variant_score", "0")),
        "structural": safe_float(row.get("structural_score", "0")),
        "dependency": safe_float(row.get("dependency_score", "0")),
        "expression": safe_float(row.get("expression_score", "0")),
        "drug": safe_float(row.get("drug_evidence_score", "0")),
    }
    score = sum(WEIGHTS[k] * components[k] for k in WEIGHTS)
    return round(score, 4), components


def render_report(
    sample_id: str,
    rows: list[dict],
    *,
    pipeline_version: str,
    input_path: Path,
    input_sha: str,
    n_variants: int,
    counts: dict[str, int],
    subtype: str,
) -> str:
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    top = rows[:10]

    lines: list[str] = []
    lines.append(f"# RMS-ISP Therapeutic Hypothesis Report")
    lines.append("")
    lines.append(f"- **Sample**: `{sample_id}`")
    lines.append(f"- **Subtype assumed**: `{subtype}`")
    lines.append(f"- **Pipeline version**: `{pipeline_version}`")
    lines.append(f"- **Run timestamp (UTC)**: {now}")
    lines.append(f"- **Input VCF**: `{input_path.name}` (sha256[:16] `{input_sha}`)")
    lines.append(f"- **Events in**: {n_variants}  ·  drivers {counts.get('DRIVER', 0)}  ·  VUS {counts.get('VUS', 0)}  ·  passengers {counts.get('PASSENGER', 0)}  ·  off-target {counts.get('OFF_TARGET', 0)}")
    lines.append("")
    lines.append("> Engineering output, not medical advice. Phase 3 dependency uses real DepMap 26Q1 Chronos with an oncogene-addiction floor for activated drivers; Phase 3 expression uses real OpenPedCan v15 RNA-seq (RMS vs other pediatric, log2 TPM z-score) when available. Phase 4 unions the curated drug map with live DGIdb interactions and ClinicalTrials.gov pediatric-RMS trial upgrades. Use to verify pipeline behaviour and triage hypotheses; do not use for clinical decisions.")
    lines.append("")
    lines.append("## Top-ranked target-drug pairs")
    lines.append("")
    lines.append("| Rank | Gene | Event | Type | Call | Drug | Mechanism | Phase | Ped | **Confidence** |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|")
    for i, r in enumerate(top, 1):
        lines.append(
            f"| {i} | **{r['gene']}** | {event_label(r)} | {r.get('event_type', '')} | {r['call']} | "
            f"`{r['drug']}` | {r['drug_mechanism']} | {r['drug_max_phase']} | "
            f"{r['drug_pediatric_evidence']} | **{r['confidence']:.3f}** |"
        )
    lines.append("")
    lines.append("## Per-event detail")
    lines.append("")
    by_event: dict[str, list[dict]] = {}
    for r in rows:
        by_event.setdefault(r["event_id"], []).append(r)
    # Sort event groups by their best confidence score (descending).
    event_groups = sorted(by_event.items(), key=lambda kv: -max(x["confidence"] for x in kv[1]))
    for vid, vrows in event_groups:
        head = vrows[0]
        lines.append(f"### `{vid}` — {head['gene']} {event_label(head)}  ({head.get('event_type', '')}, {head['call']})")
        lines.append("")
        lines.append(f"- **Reason**: {head['reason']}")
        lines.append(f"- **Variant score**: {safe_float(head['variant_score']):.2f}  ·  **Structural score**: {safe_float(head['structural_score']):.2f}  ·  **Dependency score**: {safe_float(head['dependency_score']):.2f}")
        lines.append(f"- **Dependency context**: {head.get('dependency_reason', '')}")
        if head.get("alphafold_url"):
            lines.append(f"- **Structural reference**: [{head['uniprot']}]({head['alphafold_url']})")
        lines.append("")
        lines.append("| Drug | Mechanism | Phase | Ped evidence | Drug score | **Confidence** |")
        lines.append("|---|---|---|---|---|---|")
        vrows_sorted = sorted(vrows, key=lambda r: -r["confidence"])
        for r in vrows_sorted:
            lines.append(
                f"| `{r['drug']}` | {r['drug_mechanism']} | {r['drug_max_phase']} | "
                f"{r['drug_pediatric_evidence']} | {safe_float(r['drug_evidence_score']):.2f} | "
                f"**{r['confidence']:.3f}** |"
            )
        lines.append("")

    lines.append("## Methodology")
    lines.append("")
    lines.append("Confidence formula (pilot plan §4 Aim 2 Phase 5):")
    lines.append("")
    lines.append("```")
    for k, w in WEIGHTS.items():
        lines.append(f"  {k:<12s} weight = {w:.2f}")
    lines.append("")
    lines.append("  confidence = sum(weight_k * component_k) for k in {variant, structural, dependency, expression, drug}")
    lines.append("```")
    lines.append("")
    lines.append("Component definitions:")
    lines.append("- **variant**: 1.0 for known hotspots and LoF in tumor suppressors, 0.4-0.6 for protein-altering VUS in target genes, 0 for passengers.")
    lines.append("- **structural**: 1.0 for hotspot residues with reference AlphaFold structures, 0.6 for VUS missense, 0.4 for LoF in TSGs, 0 otherwise.")
    lines.append("- **dependency**: derived from DepMap Chronos scores across RMS cell lines (mean and FP/FN-stratified), with subtype-aware bonus when the requested subtype is more dependent than the all-RMS mean.")
    lines.append("- **expression**: log2(TPM+1) z-score of RMS vs other pediatric cancer RNA-seq from OpenPedCan v15; clipped to [0,1] so a gene must be at least 2 standard deviations above the cross-disease mean to score 1.0.")
    lines.append("- **drug**: phase_weight × pediatric_evidence_weight from the curated drug-target map.")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("*Generated by `rms-isp` Phase 5. This is engineering validation output, not medical advice.*")
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, type=Path,
                    help="Phase 4 long-format TSV (one row per variant-drug pair).")
    ap.add_argument("--vcf", required=True, type=Path,
                    help="Original VCF, used for the report's input fingerprint.")
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--subtype", default="ALL")
    ap.add_argument("--pipeline-version", default=os.environ.get("RMS_ISP_VERSION", "v0.1.0-pilot"))
    ap.add_argument("--out-tsv", required=True, type=Path)
    ap.add_argument("--out-md", required=True, type=Path)
    args = ap.parse_args()

    with args.inp.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        in_cols = list(reader.fieldnames or [])

    out_cols = in_cols + ["confidence", "component_variant", "component_structural", "component_dependency", "component_expression", "component_drug"]

    for r in rows:
        score, comps = confidence(r)
        r["confidence"] = score
        for k, v in comps.items():
            r[f"component_{k}"] = f"{v:.3f}"

    rows.sort(key=lambda r: -r["confidence"])

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            r_out = dict(r)
            r_out["confidence"] = f"{r['confidence']:.4f}"
            w.writerow(r_out)

    # Re-derive variant counts and call distribution from the upstream rows
    # (the long table has duplicates per drug, so dedupe by event_id).
    seen_variants: dict[str, str] = {}
    for r in rows:
        seen_variants[r["event_id"]] = r["call"]
    counts: dict[str, int] = {}
    for call in seen_variants.values():
        counts[call] = counts.get(call, 0) + 1
    n_variants = len(seen_variants)

    input_sha = file_sha256(args.vcf)
    md = render_report(
        args.sample_id, rows,
        pipeline_version=args.pipeline_version,
        input_path=args.vcf,
        input_sha=input_sha,
        n_variants=n_variants,
        counts=counts,
        subtype=args.subtype,
    )
    args.out_md.write_text(md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
