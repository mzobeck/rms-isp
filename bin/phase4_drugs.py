#!/usr/bin/env python3
"""
Phase 4: Drug matching.

Joins each variant's gene against assets/drug_target_map.tsv and emits a
per-variant + per-drug long-format TSV. Each variant fans out into one row per
candidate drug with a drug_evidence_score.

Conditional rules honoured here (v0.1):
  - KRAS_G12C_inhibitor (adagrasib, sotorasib) is suppressed unless the
    matched variant's hgvsp_short is exactly G12C.
  - tipifarnib (HRAS-selective) is only emitted when role == "oncogene"
    and gene is HRAS, which the table already constrains.
  - MDM2 inhibitors are emitted regardless of TP53 status in v0.1; a future
    `applies_when` column will let the table express this constraint.

drug_evidence_score = phase_weight * pediatric_weight, where:
  phase_weight     = approved 1.0 | phase3 0.85 | phase2 0.70 | phase1 0.55 | preclinical 0.30
  pediatric_weight = yes_approved 1.0 | yes_trial 0.85 | adult_only 0.65 | none 0.40
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


PHASE_WEIGHT = {
    "approved": 1.0,
    "phase3": 0.85,
    "phase2": 0.70,
    "phase1": 0.55,
    "preclinical": 0.30,
}
PED_WEIGHT = {
    "yes_approved": 1.0,
    "yes_trial": 0.85,
    "adult_only": 0.65,
    "none": 0.40,
}


def load_drugs(path: Path) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    with path.open() as fh:
        reader = csv.DictReader((ln for ln in fh if not ln.startswith("##")), delimiter="\t")
        for row in reader:
            out.setdefault(row["gene"].strip(), []).append(row)
    return out


def drug_applies(drug_row: dict, variant: dict) -> bool:
    mech = drug_row.get("mechanism", "")
    hgvsp = variant.get("hgvsp_short", "")
    if mech == "KRAS_G12C_inhibitor" and hgvsp != "G12C":
        return False
    return True


def score_drug(drug_row: dict) -> float:
    p = PHASE_WEIGHT.get(drug_row.get("max_phase", ""), 0.30)
    e = PED_WEIGHT.get(drug_row.get("pediatric_evidence", ""), 0.40)
    return round(p * e, 3)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, type=Path)
    ap.add_argument("--drug-map", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    drug_map = load_drugs(args.drug_map)

    with args.inp.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        in_rows = list(reader)
        in_cols = list(reader.fieldnames or [])

    add_cols = [
        "drug", "drug_mechanism", "drug_max_phase", "drug_pediatric_evidence",
        "drug_evidence_score", "drug_notes",
    ]
    out_cols = in_cols + add_cols
    out_rows: list[dict] = []

    def emit_no_drug(v: dict, reason: str) -> dict:
        row = dict(v)
        row["drug"] = ""
        row["drug_mechanism"] = ""
        row["drug_max_phase"] = ""
        row["drug_pediatric_evidence"] = ""
        row["drug_evidence_score"] = "0.000"
        row["drug_notes"] = reason
        return row

    for v in in_rows:
        gene = v.get("gene", "")
        call = v.get("call", "")
        # Passengers and off-target variants survive into the long table with
        # an empty drug slot so downstream counts and reports stay honest.
        if call not in {"DRIVER", "VUS"}:
            out_rows.append(emit_no_drug(v, f"no drug matching: variant call is {call}"))
            continue
        candidates = drug_map.get(gene, [])
        applicable = [d for d in candidates if drug_applies(d, v)]
        if not applicable:
            out_rows.append(emit_no_drug(v, f"no drug matching: no applicable drug for {gene} in v0.1 map"))
            continue
        for drug in applicable:
            row = dict(v)
            row["drug"] = drug["drug"]
            row["drug_mechanism"] = drug["mechanism"]
            row["drug_max_phase"] = drug["max_phase"]
            row["drug_pediatric_evidence"] = drug["pediatric_evidence"]
            row["drug_evidence_score"] = f"{score_drug(drug):.3f}"
            row["drug_notes"] = drug.get("notes", "")
            out_rows.append(row)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
