#!/usr/bin/env python3
"""
Phase 2: Structural reference attachment.

For each event on a target gene with a UniProt ID, attach the AlphaFold DB
entry URL and a structural-confidence score:

  SNV (DRIVER hotspot):       1.0  reference structure resolves the residue
  SNV (DRIVER LoF in TSG):    0.4  no folded mutant; gene-level call
  SNV (VUS missense):         0.6  reference available; mutant prediction deferred
  CNA amplification (DRIVER): 0.5  gene-level event; structural model unaffected
  CNA homozygous deletion:    0.4  no protein; reference still useful for context
  Fusion (DRIVER):            0.7  fusion partner structures available; junction modelling deferred
  PASSENGER / OFF_TARGET:     0.0

v0.1 ships AlphaFold *reference* attachment only. Real Boltz-1 / Chai-1 mutant
prediction, AlphaFold-Multimer for fusion junctions, and AutoDock Vina ligand
docking are scoped for v0.3 on a shortlist.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


ALPHAFOLD_URL_TEMPLATE = "https://alphafold.ebi.ac.uk/entry/{uniprot}"
ALPHAFOLD_PDB_TEMPLATE = "https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb"


def structural_score(call: str, event_type: str, hgvsp_short: str, role: str) -> tuple[float, str]:
    if call not in {"DRIVER", "VUS"}:
        return (0.0, "no structural relevance")
    if event_type == "snv":
        if call == "DRIVER" and hgvsp_short:
            return (1.0, "known hotspot; reference structure resolves residue")
        if call == "DRIVER" and role == "tsg":
            return (0.4, "LoF in TSG; reference structure available, mutant truncated")
        return (0.6, "protein-altering missense; reference AlphaFold available; mutant prediction deferred")
    if event_type == "cna":
        if call == "DRIVER" and "deletion" in (hgvsp_short or "") or "homozygous" in role.lower():
            pass  # fall through to next check
        # Use the consequence column to disambiguate amp vs del
        # (we receive consequence in the row dict, not here — re-derive caller-side).
        return (0.5, "gene-level CNA; reference protein structure available")
    if event_type == "fusion":
        return (0.7, "fusion partner reference structures available; junction modelling deferred")
    return (0.0, "unknown event type")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    with args.inp.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        in_cols = list(reader.fieldnames or [])

    add_cols = ["alphafold_url", "alphafold_pdb_url", "structural_score", "structural_reason"]
    out_cols = in_cols + add_cols

    for r in rows:
        uniprot = r.get("uniprot", "")
        call = r.get("call", "")
        if uniprot and call != "OFF_TARGET":
            r["alphafold_url"] = ALPHAFOLD_URL_TEMPLATE.format(uniprot=uniprot)
            r["alphafold_pdb_url"] = ALPHAFOLD_PDB_TEMPLATE.format(uniprot=uniprot)
        else:
            r["alphafold_url"] = ""
            r["alphafold_pdb_url"] = ""
        score, reason = structural_score(call, r.get("event_type", ""), r.get("hgvsp_short", ""), r.get("role", ""))
        # CNA-specific refinement: deletion events are LoF-like, amplifications keep the protein.
        if r.get("event_type") == "cna" and call == "DRIVER":
            ev = (r.get("consequence") or "").lower()
            if "deletion" in ev or ev in {"homdel", "deep_deletion"}:
                score, reason = 0.4, "homozygous deletion; no protein product, reference structure for context"
            elif "amplification" in ev or ev in {"amp", "high_amp"}:
                score, reason = 0.5, "gene amplification; reference protein structure available, dosage effect"
        r["structural_score"] = f"{score:.3f}"
        r["structural_reason"] = reason

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
