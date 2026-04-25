#!/usr/bin/env python3
"""
Fetch real RMS tumor mutation calls from the TARGET-RT cohort published on
cBioPortal (study rms_nih_2014, the Shern et al. 2014 Cancer Discov cohort:
43 tumor / normal whole-genome or whole-exome pairs).

Source: https://www.cbioportal.org/api/. Public, no auth.

For each sample, writes a synthetic VCF in our toy-VCF format (with GENE,
CONSEQUENCE, NOTE in INFO) plus a per-sample row in a sample manifest.
Subtype (FP vs FN) is auto-classified from cBioPortal's PAX_FUSION clinical
attribute when available.

Outputs land under data/target_rt/ which is gitignored. Only the cohort
manifest is committed to the repo (under results/target_rt_manifest.tsv when
the runner is invoked).

Usage:
    bin/fetch_target_rt.py                    # fetch all 43 samples
    bin/fetch_target_rt.py --limit 5          # fetch first 5 with target hits

The pipeline can then be run per sample via main.nf, or all-at-once via
bin/run_target_rt.sh.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
CBIOPORTAL = "https://www.cbioportal.org/api"
STUDY_ID = "rms_nih_2014"
MUTATION_PROFILE = "rms_nih_2014_mutations"
SAMPLE_LIST = "rms_nih_2014_all"

# cBioPortal mutationType -> our pipeline's CONSEQUENCE vocabulary
CONSEQUENCE_MAP = {
    "Missense_Mutation":     "missense",
    "Nonsense_Mutation":     "stop_gained",
    "Nonstop_Mutation":      "stop_lost",
    "Frame_Shift_Del":       "frameshift",
    "Frame_Shift_Ins":       "frameshift",
    "In_Frame_Del":          "inframe_deletion",
    "In_Frame_Ins":          "inframe_insertion",
    "Splice_Site":           "splice_acceptor",
    "Silent":                "synonymous",
    "Intron":                "intron_variant",
    "5'UTR":                 "five_prime_utr_variant",
    "3'UTR":                 "three_prime_utr_variant",
    "5'Flank":               "upstream_gene_variant",
    "3'Flank":               "downstream_gene_variant",
    "IGR":                   "intergenic_variant",
    "RNA":                   "non_coding_transcript_exon_variant",
    "Translation_Start_Site":"start_lost",
}


def http_get(url: str) -> dict | list:
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read())


def http_post(url: str, body: dict) -> list:
    req = urllib.request.Request(
        url, data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read())


def load_target_entrez(targets_kb: Path, gene_effect_csv: Path) -> dict[str, str]:
    """Return {gene_symbol: entrez_id} for our 21 targets, sourced from the
    DepMap gene_effect.csv header which has 'SYMBOL (entrez_id)' columns."""
    import re
    targets: set[str] = set()
    with targets_kb.open() as fh:
        reader = csv.DictReader((ln for ln in fh if not ln.startswith("##")), delimiter="\t")
        for row in reader:
            targets.add(row["gene"].strip())
    out: dict[str, str] = {}
    with gene_effect_csv.open() as fh:
        header = fh.readline().rstrip("\n").split(",")
    for h in header[1:]:
        m = re.match(r"^(\S+)\s+\((\d+)\)$", h)
        if m and m.group(1) in targets:
            out[m.group(1)] = m.group(2)
    return out


def fetch_samples(study: str) -> list[dict]:
    return http_get(f"{CBIOPORTAL}/studies/{study}/samples")


def fetch_clinical(study: str) -> dict[str, dict]:
    """Return {sample_id: {attribute_id: value}}."""
    rows = http_get(f"{CBIOPORTAL}/studies/{study}/clinical-data?clinicalDataType=SAMPLE")
    out: dict[str, dict] = {}
    for r in rows:
        sid = r.get("sampleId", "")
        out.setdefault(sid, {})[r["clinicalAttributeId"]] = r["value"]
    return out


def fetch_mutations(profile: str, sample_list: str, entrez_ids: list[str]) -> list[dict]:
    return http_post(
        f"{CBIOPORTAL}/molecular-profiles/{profile}/mutations/fetch?projection=DETAILED",
        {"sampleListId": sample_list, "entrezGeneIds": [int(e) for e in entrez_ids]},
    )


def classify_subtype(clinical: dict) -> str:
    """Map cBioPortal PAX_FUSION value to our --subtype param."""
    pax = (clinical.get("PAX_FUSION") or "").upper()
    if "PAX3" in pax or "PAX7" in pax or "POSITIVE" in pax or pax == "FP":
        return "FP"
    if "NEGATIVE" in pax or "NONE" in pax or pax == "FN":
        return "FN"
    return "ALL"


def write_vcf(sample_id: str, mutations: list[dict], out_path: Path) -> int:
    """Write a per-sample VCF in our toy-VCF format. Returns mutation count."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    lines.append("##fileformat=VCFv4.2")
    lines.append(f"##fileDate={datetime.now(timezone.utc).strftime('%Y%m%d')}")
    lines.append(f"##source=RMS-ISP fetch_target_rt.py from cBioPortal study {STUDY_ID}")
    lines.append("##reference=GRCh37 (cBioPortal/TARGET-RT default; coordinates emitted as-is)")
    lines.append('##INFO=<ID=GENE,Number=1,Type=String,Description="Affected gene symbol">')
    lines.append('##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description="Variant consequence">')
    lines.append('##INFO=<ID=NOTE,Number=1,Type=String,Description="Curator note (carries proteinChange for hotspot extraction)">')
    lines.append('##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic variant">')
    lines.append('##INFO=<ID=SOURCE,Number=1,Type=String,Description="Original cBioPortal mutation row source">')
    lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    lines.append(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}")
    for i, m in enumerate(mutations, 1):
        chrom = f"chr{m.get('chr','')}" if not str(m.get('chr','')).startswith("chr") else m.get('chr','')
        pos = m.get("startPosition", "")
        ref = m.get("referenceAllele") or "."
        alt = m.get("variantAllele") or "."
        gene = (m.get("gene") or {}).get("hugoGeneSymbol", "")
        mtype = m.get("mutationType", "")
        consequence = CONSEQUENCE_MAP.get(mtype, mtype.lower() or "unknown")
        protein_change = m.get("proteinChange") or ""
        # Strip leading "p." if present so Phase 1's regex picks the AA change cleanly
        protein_change = protein_change.lstrip("p.").strip()
        note = f"{protein_change}_targetrt_{m.get('mutationStatus','SOMATIC')}".replace(" ", "_") if protein_change else f"targetrt_{mtype}"
        info = f"GENE={gene};CONSEQUENCE={consequence};NOTE={note};SOMATIC;SOURCE=cbioportal"
        lines.append(f"{chrom}\t{pos}\trms2014_var{i:03d}\t{ref}\t{alt}\t.\tPASS\t{info}\tGT\t0/1")
    out_path.write_text("\n".join(lines) + "\n")
    return len(mutations)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", type=Path,
                    default=REPO_ROOT / "data" / "target_rt",
                    help="Per-sample VCF output directory (gitignored).")
    ap.add_argument("--manifest", type=Path,
                    default=REPO_ROOT / "data" / "target_rt" / "manifest.tsv")
    ap.add_argument("--targets-kb", type=Path,
                    default=REPO_ROOT / "assets" / "targets_kb.tsv")
    ap.add_argument("--gene-effect", type=Path, default=Path("/tmp/depmap/gene_effect.csv"),
                    help="Used for the SYMBOL -> Entrez ID mapping.")
    ap.add_argument("--limit", type=int, default=0,
                    help="If >0, only emit the first N samples that have at least one target hit.")
    ap.add_argument("--include-untargeted", action="store_true",
                    help="Also emit VCFs for samples with zero hits on our 21 targets.")
    args = ap.parse_args()

    # 1. Get Entrez IDs for our targets (from DepMap header).
    entrez = load_target_entrez(args.targets_kb, args.gene_effect)
    print(f"loaded Entrez IDs for {len(entrez)} target genes", file=sys.stderr)

    # 2. Pull samples + clinical + mutations.
    samples = fetch_samples(STUDY_ID)
    print(f"study {STUDY_ID}: {len(samples)} samples", file=sys.stderr)
    clinical = fetch_clinical(STUDY_ID)
    print(f"clinical attributes loaded for {len(clinical)} samples", file=sys.stderr)
    mutations = fetch_mutations(MUTATION_PROFILE, SAMPLE_LIST, list(entrez.values()))
    print(f"fetched {len(mutations)} mutation rows hitting our 21 targets", file=sys.stderr)

    # 3. Bucket mutations by sample.
    by_sample: dict[str, list[dict]] = {}
    for m in mutations:
        by_sample.setdefault(m["sampleId"], []).append(m)

    # 4. Emit VCFs + manifest.
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)

    sample_order = sorted(samples, key=lambda s: -len(by_sample.get(s["sampleId"], [])))
    if args.limit > 0:
        sample_order = [s for s in sample_order if by_sample.get(s["sampleId"])][:args.limit]
    elif not args.include_untargeted:
        sample_order = [s for s in sample_order if by_sample.get(s["sampleId"])]

    manifest_rows: list[dict] = []
    for s in sample_order:
        sid = s["sampleId"]
        muts = by_sample.get(sid, [])
        clin = clinical.get(sid, {})
        subtype = classify_subtype(clin)
        vcf_path = args.out_dir / f"{sid}.vcf"
        n = write_vcf(sid, muts, vcf_path)
        manifest_rows.append({
            "sample_id": sid,
            "subtype": subtype,
            "pax_fusion": clin.get("PAX_FUSION", ""),
            "histology": clin.get("HISTOLOGICAL_SUBTYPE", ""),
            "risk_group": clin.get("RISK_GROUP", ""),
            "tumor_location": clin.get("PRIMARY_TUMOR_LOCATION", ""),
            "age": clin.get("AGE", ""),
            "sex": clin.get("SEX", ""),
            "tmb": clin.get("TMB_NONSYNONYMOUS", ""),
            "n_mutations_on_targets": n,
            "vcf_path": str(vcf_path.relative_to(REPO_ROOT)),
        })

    cols = list(manifest_rows[0].keys()) if manifest_rows else []
    with args.manifest.open("w", newline="") as fh:
        fh.write(f"## TARGET-RT cohort manifest (auto-generated by bin/fetch_target_rt.py)\n")
        fh.write(f"## Source: cBioPortal study {STUDY_ID} (Shern 2014 Cancer Discov)\n")
        fh.write(f"## Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}\n")
        fh.write(f"## Samples written: {len(manifest_rows)}\n")
        if cols:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(manifest_rows)
    print(f"wrote {len(manifest_rows)} VCFs + {args.manifest}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
