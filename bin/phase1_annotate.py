#!/usr/bin/env python3
"""
Phase 1: Variant annotation.

Reads a VCF, joins each variant against assets/targets_kb.tsv to determine
target-list membership and driver/passenger status, and writes a per-variant
TSV consumed by downstream phases.

v0.1 leans on the GENE/CONSEQUENCE INFO fields supplied in the toy VCF. A real
Phase 1 will run VEP and OncoKB; the column contract emitted here matches what
that downstream module would produce.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


HOTSPOT_RE = re.compile(r"(?<![A-Za-z0-9])([A-Z]\d+[A-Z*X])(?![A-Za-z0-9])")
PROTEIN_ALTERING_CONSEQUENCES = {
    "missense",
    "missense_variant",
    "stop_gained",
    "stop_lost",
    "frameshift",
    "frameshift_variant",
    "inframe_insertion",
    "inframe_deletion",
    "splice_acceptor",
    "splice_donor",
    "splice_acceptor_variant",
    "splice_donor_variant",
}
LOSS_OF_FUNCTION_CONSEQUENCES = {
    "stop_gained",
    "frameshift",
    "frameshift_variant",
    "splice_acceptor",
    "splice_donor",
    "splice_acceptor_variant",
    "splice_donor_variant",
}


def parse_info(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for token in info.split(";"):
        if "=" in token:
            k, v = token.split("=", 1)
            out[k] = v
        else:
            out[token] = "1"
    return out


def load_targets_kb(path: Path) -> dict[str, dict]:
    kb: dict[str, dict] = {}
    with path.open() as fh:
        reader = csv.DictReader((ln for ln in fh if not ln.startswith("##")), delimiter="\t")
        for row in reader:
            gene = row["gene"].strip()
            kb[gene] = {
                "uniprot": row["uniprot"].strip(),
                "role": row["role"].strip(),
                "hotspots": {h.strip() for h in row["hotspots"].split(",") if h.strip()},
                "lof_target": row["loss_of_function_target"].strip() == "1",
                "notes": row["notes"].strip(),
            }
    return kb


def extract_hotspot(note: str) -> str:
    m = HOTSPOT_RE.search(note or "")
    return m.group(1) if m else ""


def classify(consequence: str, gene_kb: dict | None, hgvsp_short: str) -> tuple[str, str, float]:
    """Return (call, reason, variant_score in 0..1)."""
    if gene_kb is None:
        return ("OFF_TARGET", "gene not in target KB", 0.0)
    cons = (consequence or "").lower()
    is_protein_altering = cons in PROTEIN_ALTERING_CONSEQUENCES
    is_lof = cons in LOSS_OF_FUNCTION_CONSEQUENCES
    if is_lof and gene_kb["lof_target"]:
        return ("DRIVER", f"LoF in TSG ({cons})", 1.0)
    if is_protein_altering and hgvsp_short and hgvsp_short in gene_kb["hotspots"]:
        return ("DRIVER", f"hotspot {hgvsp_short}", 1.0)
    if is_protein_altering and gene_kb["hotspots"]:
        return ("VUS", f"protein-altering ({cons}) in target gene; not a known hotspot", 0.4)
    if is_protein_altering:
        return ("VUS", f"protein-altering ({cons}) in target gene; no hotspot list", 0.5)
    return ("PASSENGER", f"non-protein-altering consequence ({cons})", 0.0)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, type=Path)
    ap.add_argument("--targets-kb", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--sample-id", default="TOY_TUMOR")
    args = ap.parse_args()

    kb = load_targets_kb(args.targets_kb)

    cols = [
        "sample_id", "variant_id", "chrom", "pos", "ref", "alt",
        "gene", "uniprot", "role", "consequence", "hgvsp_short",
        "is_target", "call", "reason", "variant_score",
    ]
    rows: list[dict[str, str]] = []
    with args.vcf.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, vid, ref, alt, _qual, _filt, info = parts[:8]
            kv = parse_info(info)
            gene = kv.get("GENE", "")
            consequence = kv.get("CONSEQUENCE", "")
            note = kv.get("NOTE", "")
            hgvsp_short = extract_hotspot(note)
            gene_kb = kb.get(gene)
            call, reason, vscore = classify(consequence, gene_kb, hgvsp_short)
            rows.append({
                "sample_id": args.sample_id,
                "variant_id": vid,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "gene": gene,
                "uniprot": gene_kb["uniprot"] if gene_kb else "",
                "role": gene_kb["role"] if gene_kb else "",
                "consequence": consequence,
                "hgvsp_short": hgvsp_short,
                "is_target": "1" if gene_kb else "0",
                "call": call,
                "reason": reason,
                "variant_score": f"{vscore:.3f}",
            })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    n_driver = sum(1 for r in rows if r["call"] == "DRIVER")
    n_vus = sum(1 for r in rows if r["call"] == "VUS")
    n_passenger = sum(1 for r in rows if r["call"] == "PASSENGER")
    n_off = sum(1 for r in rows if r["call"] == "OFF_TARGET")
    print(
        f"phase1: {len(rows)} variants -> {n_driver} DRIVER, {n_vus} VUS, "
        f"{n_passenger} PASSENGER, {n_off} OFF_TARGET",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
