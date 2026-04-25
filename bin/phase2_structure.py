#!/usr/bin/env python3
"""
Phase 2: Structural reference attachment.

For each Phase 1 variant on a target gene with a UniProt ID, attach the
AlphaFold DB entry URL and a structural-confidence score:
  - 1.0 if the variant is a known hotspot (binding-relevant by curation)
  - 0.6 if protein-altering missense at a non-hotspot residue (potential pocket impact)
  - 0.4 if a LoF event in a TSG (no folded mutant protein, but loss is interpretable)
  - 0.0 otherwise (passenger / off-target)

v0.1 ships AlphaFold *reference* attachment only. Real Boltz-1 / Chai-1 mutant
prediction and AutoDock Vina ligand docking are scoped for v0.2 on a shortlist.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


ALPHAFOLD_URL_TEMPLATE = "https://alphafold.ebi.ac.uk/entry/{uniprot}"
ALPHAFOLD_PDB_TEMPLATE = "https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb"


def structural_score(call: str, hgvsp_short: str, role: str) -> tuple[float, str]:
    if call == "DRIVER" and hgvsp_short:
        return (1.0, "known hotspot; reference structure available")
    if call == "DRIVER" and role == "tsg":
        return (0.4, "LoF in TSG; structural model available for reference but mutant protein is truncated")
    if call == "VUS":
        return (0.6, "protein-altering missense in target; reference AlphaFold available; mutant prediction deferred")
    return (0.0, "no structural relevance")


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
        hgvsp = r.get("hgvsp_short", "")
        role = r.get("role", "")
        if uniprot and call != "OFF_TARGET":
            r["alphafold_url"] = ALPHAFOLD_URL_TEMPLATE.format(uniprot=uniprot)
            r["alphafold_pdb_url"] = ALPHAFOLD_PDB_TEMPLATE.format(uniprot=uniprot)
        else:
            r["alphafold_url"] = ""
            r["alphafold_pdb_url"] = ""
        score, reason = structural_score(call, hgvsp, role)
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
