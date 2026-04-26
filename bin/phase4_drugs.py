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


def union_drug_maps(primary: Path, extras: list[Path]) -> dict[str, list[dict]]:
    """Union the primary curated drug map with optional extra sources (e.g.,
    DGIdb cache). Primary always wins on (gene, drug.lower()) collisions
    because its mechanism strings drive the case-study scorecard assertions.
    """
    merged = load_drugs(primary)
    seen: set[tuple[str, str]] = {
        (g, (r.get("drug") or "").strip().lower())
        for g, rows in merged.items() for r in rows
    }
    for path in extras:
        if not path or not path.exists() or path.stat().st_size == 0:
            continue
        for gene, rows in load_drugs(path).items():
            for r in rows:
                key = (gene, (r.get("drug") or "").strip().lower())
                if key in seen:
                    continue
                seen.add(key)
                merged.setdefault(gene, []).append(r)
    return merged


def load_ctgov_trials(path: Path | None) -> list[dict]:
    """Load the CT.gov pediatric-RMS trial map (output of bin/fetch_clinicaltrials.py).

    Returns a list of dicts; matching is done via substring later because
    CT.gov drug names sometimes carry suffixes ('selumetinib sulfate') or
    prefixes ('apo-trametinib') that exact match misses.
    """
    if not path or not path.exists() or path.stat().st_size == 0:
        return []
    out: list[dict] = []
    with path.open() as fh:
        reader = csv.DictReader((ln for ln in fh if not ln.startswith("##")), delimiter="\t")
        for row in reader:
            if row.get("any_pediatric", "0") != "1":
                continue
            out.append({
                "drug_norm": (row.get("drug") or "").strip().lower(),
                "n_trials": row.get("n_rms_trials", ""),
                "max_phase": row.get("max_phase_in_rms", ""),
                "any_recruiting": row.get("any_recruiting", "0") == "1",
                "example_nct": row.get("example_nct", ""),
            })
    return out


def ctgov_lookup(drug_name: str, trials: list[dict]) -> dict | None:
    """Find a CT.gov entry whose drug name overlaps with the given drug.

    We match by either substring direction so generic-name variants
    (e.g., 'selumetinib' vs 'selumetinib sulfate') still hit.
    """
    n = drug_name.strip().lower()
    if not n:
        return None
    for t in trials:
        ct_name = t["drug_norm"]
        if not ct_name:
            continue
        if ct_name == n or ct_name in n or n in ct_name:
            return t
    return None


def apply_ctgov_upgrade(row: dict, trials: list[dict]) -> None:
    """In-place upgrade of pediatric_evidence and drug_notes based on CT.gov.

    Rule: drug found in pediatric RMS trial -> pediatric_evidence becomes
    'yes_trial', UNLESS already 'yes_approved' (which outranks).
    The NCT ID is appended to drug_notes for auditability.
    """
    drug = row.get("drug", "")
    if not drug:
        return
    hit = ctgov_lookup(drug, trials)
    if not hit:
        return
    current = row.get("drug_pediatric_evidence", "")
    if current != "yes_approved":
        row["drug_pediatric_evidence"] = "yes_trial"
    note = f"ctgov:{hit['example_nct']} ({hit['n_trials']} RMS trials"
    if hit["any_recruiting"]:
        note += ", recruiting"
    note += ")"
    existing = row.get("drug_notes", "")
    row["drug_notes"] = f"{existing} | {note}" if existing else note


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
    ap.add_argument("--drug-map", required=True, type=Path,
                    help="Primary curated drug-target map.")
    ap.add_argument("--drug-map-extra", action="append", default=[], type=Path,
                    help="Optional extra drug-target maps (e.g. DGIdb cache); curated wins on (gene, drug) collisions.")
    ap.add_argument("--ctgov-trials", type=Path, default=None,
                    help="Optional CT.gov pediatric RMS trial map (output of bin/fetch_clinicaltrials.py); upgrades pediatric_evidence for drugs found in pediatric RMS trials.")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    drug_map = union_drug_maps(args.drug_map, args.drug_map_extra)
    ctgov_trials = load_ctgov_trials(args.ctgov_trials)

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
            row["drug_notes"] = drug.get("notes", "")
            # CT.gov upgrade can promote pediatric_evidence; do it BEFORE scoring
            # so the upgraded value flows into drug_evidence_score.
            apply_ctgov_upgrade(row, ctgov_trials)
            row["drug_evidence_score"] = f"{score_drug({**drug, 'pediatric_evidence': row['drug_pediatric_evidence']}):.3f}"
            out_rows.append(row)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
