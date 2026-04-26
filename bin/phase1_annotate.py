#!/usr/bin/env python3
"""
Phase 1: Variant + CNA + fusion annotation.

Reads up to three input files (VCF, CNA TSV, fusion TSV) and emits a single
unified per-event TSV with a stable column contract for downstream phases.

Each event is classified DRIVER / VUS / PASSENGER / OFF_TARGET by joining
against assets/targets_kb.tsv.

v0.1 leaned on the GENE/CONSEQUENCE INFO fields from the toy VCF; v0.2 keeps
that for SNVs and adds two additional input formats:

  CNA TSV:    sample_id\tgene\tevent\tcopy_number\tnotes
              event ∈ {amplification, homozygous_deletion, focal_gain, focal_loss}
  Fusion TSV: sample_id\tgene_5p\tgene_3p\tfusion_name\tfusion_class\tnotes

Real Phase 1 will run VEP + OncoKB for SNVs, GISTIC2 / cnvkit / FACETS for CNAs,
and STAR-Fusion / Arriba for fusions. The output column contract here is what
that real implementation emits.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

# Make sibling annotators package importable when run as a script.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from annotators import Variant, get_annotator  # noqa: E402
from annotators.curated_vcf import HOTSPOT_RE, extract_hotspot  # noqa: E402  # re-exported for backwards compat
PROTEIN_ALTERING_CONSEQUENCES = {
    "missense", "missense_variant",
    "stop_gained", "stop_lost",
    "frameshift", "frameshift_variant",
    "inframe_insertion", "inframe_deletion",
    "splice_acceptor", "splice_donor",
    "splice_acceptor_variant", "splice_donor_variant",
}
LOSS_OF_FUNCTION_CONSEQUENCES = {
    "stop_gained",
    "frameshift", "frameshift_variant",
    "splice_acceptor", "splice_donor",
    "splice_acceptor_variant", "splice_donor_variant",
}

OUT_COLS = [
    "sample_id", "event_id", "event_type", "chrom", "pos", "ref", "alt",
    "gene", "uniprot", "role",
    "consequence", "hgvsp_short", "copy_number", "fusion_partner",
    "is_target", "call", "reason", "variant_score",
]


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
                "amp_target": row.get("oncogenic_amplification", "0").strip() == "1",
                "del_target": row.get("lof_via_deletion", "0").strip() == "1",
                "notes": row["notes"].strip(),
            }
    return kb


def classify_snv(consequence: str, gene_kb: dict | None, hgvsp_short: str) -> tuple[str, str, float]:
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


def classify_cna(event: str, gene_kb: dict | None) -> tuple[str, str, float]:
    if gene_kb is None:
        return ("OFF_TARGET", "gene not in target KB", 0.0)
    ev = (event or "").lower()
    if ev in {"amplification", "amp", "high_amp"} and gene_kb["amp_target"]:
        return ("DRIVER", f"oncogenic amplification of {gene_kb['role']}", 1.0)
    if ev in {"homozygous_deletion", "homdel", "deep_deletion"} and gene_kb["del_target"]:
        return ("DRIVER", f"homozygous deletion of TSG/dependency-creating gene", 1.0)
    if ev in {"focal_gain", "gain"} and gene_kb["amp_target"]:
        return ("VUS", f"low-level gain of amp-target gene; subthreshold for DRIVER", 0.4)
    if ev in {"focal_loss", "loss"} and gene_kb["del_target"]:
        return ("VUS", f"heterozygous loss of TSG; haploinsufficiency uncertain", 0.4)
    return ("PASSENGER", f"CNA event '{ev}' on {gene_kb['role']} gene; no driver criterion met", 0.0)


def classify_fusion(gene_role: str | None, partner: str, gene_kb: dict | None) -> tuple[str, str, float]:
    if gene_kb is None:
        return ("OFF_TARGET", "gene not in target KB", 0.0)
    if gene_kb["role"] == "fusion_partner":
        return ("DRIVER", f"driver fusion partner with {partner}", 1.0)
    return ("VUS", f"fusion involving {gene_kb['role']} gene; biological role uncertain", 0.4)


def annotate_vcf(vcf_path: Path, sample_id: str, kb: dict[str, dict],
                 annotator) -> list[dict]:
    """Parse a VCF, run it through the annotator, classify each SNV.

    The annotator fills in (gene, consequence, hgvsp_short); classify_snv()
    then turns those into a DRIVER/VUS/PASSENGER call. Output column contract
    is unchanged from v0.11.
    """
    parsed: list[tuple[Variant, str]] = []  # (variant, vid)
    with vcf_path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, vid, ref, alt, _qual, _filt, info = parts[:8]
            try:
                pos_int = int(pos)
            except ValueError:
                continue
            parsed.append((
                Variant(chrom=chrom, pos=pos_int, ref=ref, alt=alt,
                        info=parse_info(info)),
                vid,
            ))

    annotations = annotator.annotate_batch([v for v, _ in parsed])

    rows: list[dict] = []
    for (v, vid), ann in zip(parsed, annotations):
        gene = ann.gene
        gene_kb = kb.get(gene)
        call, reason, vscore = classify_snv(ann.consequence, gene_kb, ann.hgvsp_short)
        rows.append({
            "sample_id": sample_id, "event_id": vid, "event_type": "snv",
            "chrom": v.chrom, "pos": str(v.pos), "ref": v.ref, "alt": v.alt,
            "gene": gene,
            "uniprot": gene_kb["uniprot"] if gene_kb else "",
            "role": gene_kb["role"] if gene_kb else "",
            "consequence": ann.consequence, "hgvsp_short": ann.hgvsp_short,
            "copy_number": "", "fusion_partner": "",
            "is_target": "1" if gene_kb else "0",
            "call": call, "reason": reason,
            "variant_score": f"{vscore:.3f}",
        })
    return rows


def annotate_cna(cna_path: Path, sample_id: str, kb: dict[str, dict]) -> list[dict]:
    rows: list[dict] = []
    with cna_path.open() as fh:
        reader = csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")
        for i, r in enumerate(reader, 1):
            if r.get("sample_id") and r["sample_id"] != sample_id:
                continue
            gene = r.get("gene", "").strip()
            event = r.get("event", "").strip()
            cn = r.get("copy_number", "").strip()
            gene_kb = kb.get(gene)
            call, reason, vscore = classify_cna(event, gene_kb)
            rows.append({
                "sample_id": sample_id, "event_id": f"cna_{i:03d}", "event_type": "cna",
                "chrom": "", "pos": "", "ref": "", "alt": "",
                "gene": gene,
                "uniprot": gene_kb["uniprot"] if gene_kb else "",
                "role": gene_kb["role"] if gene_kb else "",
                "consequence": event, "hgvsp_short": "",
                "copy_number": cn, "fusion_partner": "",
                "is_target": "1" if gene_kb else "0",
                "call": call, "reason": reason,
                "variant_score": f"{vscore:.3f}",
            })
    return rows


def annotate_fusion(fusion_path: Path, sample_id: str, kb: dict[str, dict]) -> list[dict]:
    rows: list[dict] = []
    with fusion_path.open() as fh:
        reader = csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")
        for i, r in enumerate(reader, 1):
            if r.get("sample_id") and r["sample_id"] != sample_id:
                continue
            g5 = r.get("gene_5p", "").strip()
            g3 = r.get("gene_3p", "").strip()
            name = r.get("fusion_name", "").strip() or f"{g5}-{g3}"
            for gene, partner in ((g5, g3), (g3, g5)):
                if not gene:
                    continue
                gene_kb = kb.get(gene)
                call, reason, vscore = classify_fusion(
                    gene_kb["role"] if gene_kb else None, partner, gene_kb,
                )
                rows.append({
                    "sample_id": sample_id, "event_id": f"fusion_{i:03d}_{gene}",
                    "event_type": "fusion",
                    "chrom": "", "pos": "", "ref": "", "alt": "",
                    "gene": gene,
                    "uniprot": gene_kb["uniprot"] if gene_kb else "",
                    "role": gene_kb["role"] if gene_kb else "",
                    "consequence": "fusion", "hgvsp_short": "",
                    "copy_number": "", "fusion_partner": partner,
                    "is_target": "1" if gene_kb else "0",
                    "call": call, "reason": f"{reason} ({name})",
                    "variant_score": f"{vscore:.3f}",
                })
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, type=Path)
    ap.add_argument("--cna", type=Path, default=None,
                    help="Optional CNA TSV; pass an empty file or omit if no CNAs.")
    ap.add_argument("--fusion", type=Path, default=None,
                    help="Optional fusion TSV; pass an empty file or omit if no fusions.")
    ap.add_argument("--targets-kb", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--sample-id", default="TOY_TUMOR")
    ap.add_argument("--annotator", choices=["curated", "vep_rest"], default="curated",
                    help="Variant annotator backend. 'curated' reads gene/"
                         "consequence from VCF INFO (toys, cBioPortal output). "
                         "'vep_rest' calls Ensembl VEP REST.")
    ap.add_argument("--vep-cache-dir", type=Path, default=None,
                    help="Cache dir for vep_rest annotator. "
                         "Default: <repo>/data/vep_cache.")
    ap.add_argument("--disease", default="RMS",
                    help="Reserved for future OncoKB / AlphaMissense backends; "
                         "currently ignored.")
    args = ap.parse_args()

    annotator = get_annotator(
        args.annotator,
        cache_dir=args.vep_cache_dir,
        disease=args.disease,
    )

    kb = load_targets_kb(args.targets_kb)
    rows: list[dict] = []
    rows.extend(annotate_vcf(args.vcf, args.sample_id, kb, annotator))
    if args.cna and args.cna.exists() and args.cna.stat().st_size > 0:
        rows.extend(annotate_cna(args.cna, args.sample_id, kb))
    if args.fusion and args.fusion.exists() and args.fusion.stat().st_size > 0:
        rows.extend(annotate_fusion(args.fusion, args.sample_id, kb))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    counts: dict[str, int] = {}
    for r in rows:
        counts[r["call"]] = counts.get(r["call"], 0) + 1
    by_type: dict[str, int] = {}
    for r in rows:
        by_type[r["event_type"]] = by_type.get(r["event_type"], 0) + 1
    print(
        f"phase1: {len(rows)} events ({', '.join(f'{k}={v}' for k, v in by_type.items())}) "
        f"-> {', '.join(f'{k}={v}' for k, v in sorted(counts.items()))}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
