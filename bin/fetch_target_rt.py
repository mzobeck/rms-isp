#!/usr/bin/env python3
"""
Fetch real RMS tumor data from cBioPortal and write per-sample inputs in our
pipeline's VCF + CNA TSV + fusion TSV format.

Default cohort: rms_nih_2014 (Shern 2014 Cancer Discov, 43 tumor/normal WGS/WES
pairs) plus rms_msk_2023 (MSK-IMPACT targeted sequencing of 24 extremity RMS
cases with CNA and structural-variant calls).

For each sample, writes whichever of {VCF, CNA TSV, fusion TSV} the source
study has. Subtype is auto-classified from the PAX_FUSION clinical attribute
when present, otherwise from any PAX-FOXO1 fusion call (FP) vs absence of
fusion calls (defaults to ALL).

Outputs land under data/target_rt/ which is gitignored. The combined cohort
manifest is what bin/run_target_rt.py reads.

Usage:
    bin/fetch_target_rt.py
    bin/fetch_target_rt.py --studies rms_nih_2014,rms_msk_2023
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
CBIOPORTAL = "https://www.cbioportal.org/api"

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

# cBioPortal CNA discrete alteration codes -> our event vocabulary
CNA_MAP = {
    -2: "homozygous_deletion",
    -1: "focal_loss",
    1:  "focal_gain",
    2:  "amplification",
}


def http_get(url: str) -> dict | list:
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read())


def http_post(url: str, body: dict) -> list | dict:
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


def list_profiles(study: str) -> dict[str, str]:
    """Return {alteration_type: profile_id} for the study."""
    out: dict[str, str] = {}
    for p in http_get(f"{CBIOPORTAL}/studies/{study}/molecular-profiles"):
        out[p.get("molecularAlterationType", "")] = p.get("molecularProfileId", "")
    return out


def fetch_samples(study: str) -> list[dict]:
    return http_get(f"{CBIOPORTAL}/studies/{study}/samples")


def fetch_clinical(study: str) -> dict[str, dict]:
    rows = http_get(f"{CBIOPORTAL}/studies/{study}/clinical-data?clinicalDataType=SAMPLE")
    out: dict[str, dict] = {}
    for r in rows:
        sid = r.get("sampleId", "")
        out.setdefault(sid, {})[r["clinicalAttributeId"]] = r["value"]
    return out


def fetch_mutations(profile: str, study: str, entrez_ids: list[str]) -> list[dict]:
    return http_post(
        f"{CBIOPORTAL}/molecular-profiles/{profile}/mutations/fetch?projection=DETAILED",
        {"sampleListId": f"{study}_all", "entrezGeneIds": [int(e) for e in entrez_ids]},
    )


def fetch_cnas(profile: str, study: str, entrez_ids: list[str]) -> list[dict]:
    return http_post(
        f"{CBIOPORTAL}/molecular-profiles/{profile}/discrete-copy-number/fetch?discreteCopyNumberEventType=ALL&projection=DETAILED",
        {"sampleListId": f"{study}_all", "entrezGeneIds": [int(e) for e in entrez_ids]},
    )


def fetch_svs(profile: str) -> list[dict]:
    return http_post(
        f"{CBIOPORTAL}/structural-variant/fetch",
        {"molecularProfileIds": [profile]},
    )


def classify_subtype(clinical: dict, fusion_calls: list[dict]) -> str:
    """First check PAX_FUSION clinical attribute, then fall back to fusion calls."""
    pax = (clinical.get("PAX_FUSION") or "").upper()
    if "PAX3" in pax or "PAX7" in pax or "POSITIVE" in pax or pax == "FP":
        return "FP"
    if "NEGATIVE" in pax or pax == "FN":
        return "FN"
    # Fall back to fusion calls: any PAX-FOXO1 -> FP
    for sv in fusion_calls:
        s1 = (sv.get("site1HugoSymbol") or "").upper()
        s2 = (sv.get("site2HugoSymbol") or "").upper()
        if {s1, s2} & {"PAX3", "PAX7"} and ("FOXO1" in {s1, s2}):
            return "FP"
    return "ALL"


def write_vcf(sample_id: str, mutations: list[dict], out_path: Path) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = [
        "##fileformat=VCFv4.2",
        f"##fileDate={datetime.now(timezone.utc).strftime('%Y%m%d')}",
        f"##source=RMS-ISP fetch_target_rt.py from cBioPortal",
        "##reference=GRCh37 (cBioPortal default; coordinates emitted as-is)",
        '##INFO=<ID=GENE,Number=1,Type=String,Description="Affected gene symbol">',
        '##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description="Variant consequence">',
        '##INFO=<ID=NOTE,Number=1,Type=String,Description="Curator note (carries proteinChange)">',
        '##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic variant">',
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}",
    ]
    for i, m in enumerate(mutations, 1):
        chrom = f"chr{m.get('chr','')}" if not str(m.get('chr','')).startswith("chr") else m.get('chr','')
        pos = m.get("startPosition", "")
        ref = m.get("referenceAllele") or "."
        alt = m.get("variantAllele") or "."
        gene = (m.get("gene") or {}).get("hugoGeneSymbol", "")
        mtype = m.get("mutationType", "")
        consequence = CONSEQUENCE_MAP.get(mtype, mtype.lower() or "unknown")
        protein_change = (m.get("proteinChange") or "").lstrip("p.").strip()
        note = (f"{protein_change}_targetrt" if protein_change else f"targetrt_{mtype}").replace(" ", "_")
        info = f"GENE={gene};CONSEQUENCE={consequence};NOTE={note};SOMATIC"
        lines.append(f"{chrom}\t{pos}\trms_var{i:03d}\t{ref}\t{alt}\t.\tPASS\t{info}\tGT\t0/1")
    out_path.write_text("\n".join(lines) + "\n")
    return len(mutations)


def write_cna_tsv(sample_id: str, cnas: list[dict], out_path: Path) -> int:
    """Write CNA TSV in the format Phase 1 expects.

    Filters to events where alteration is in CNA_MAP (i.e., -2/-1/+1/+2; drops 0).
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["## RMS-ISP CNA calls fetched from cBioPortal",
             "sample_id\tgene\tevent\tcopy_number\tnotes"]
    n = 0
    for c in cnas:
        alt = c.get("alteration")
        if alt not in CNA_MAP:
            continue
        gene = (c.get("gene") or {}).get("hugoGeneSymbol", "")
        if not gene:
            continue
        event = CNA_MAP[alt]
        # cBioPortal discrete CN doesn't carry the absolute copy number
        # for amplifications; record the alteration code as a proxy.
        cn_proxy = {-2: "0", -1: "1", 1: "3", 2: "5"}.get(alt, "")
        lines.append(f"{sample_id}\t{gene}\t{event}\t{cn_proxy}\tcbioportal_alt={alt}")
        n += 1
    out_path.write_text("\n".join(lines) + "\n")
    return n


def write_fusion_tsv(sample_id: str, svs: list[dict], out_path: Path) -> int:
    """Write fusion TSV in the format Phase 1 expects."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["## RMS-ISP fusion calls fetched from cBioPortal",
             "sample_id\tgene_5p\tgene_3p\tfusion_name\tfusion_class\tnotes"]
    n = 0
    for sv in svs:
        g5 = sv.get("site1HugoSymbol", "")
        g3 = sv.get("site2HugoSymbol", "")
        if not g5 or not g3 or g5 == g3:
            continue  # skip intragenic SVs (e.g. MET-MET)
        name = f"{g5}-{g3}"
        cls = sv.get("variantClass") or "fusion"
        notes = (sv.get("eventInfo") or "").replace("\n", " ")[:120]
        lines.append(f"{sample_id}\t{g5}\t{g3}\t{name}\t{cls}\t{notes}")
        n += 1
    out_path.write_text("\n".join(lines) + "\n")
    return n


def fetch_one_study(study: str, entrez: dict[str, str], out_dir: Path) -> list[dict]:
    """Fetch all available data types from one cBioPortal study and emit
    per-sample VCF + CNA TSV + fusion TSV. Returns manifest rows."""
    print(f"\n=== {study} ===", file=sys.stderr)
    profiles = list_profiles(study)
    has_mut = "MUTATION_EXTENDED" in profiles
    has_cna = "COPY_NUMBER_ALTERATION" in profiles
    has_sv  = "STRUCTURAL_VARIANT" in profiles
    print(f"  profiles: mut={has_mut} cna={has_cna} sv={has_sv}", file=sys.stderr)

    samples = fetch_samples(study)
    clinical = fetch_clinical(study)

    muts: list[dict] = []
    if has_mut:
        muts = fetch_mutations(profiles["MUTATION_EXTENDED"], study, list(entrez.values()))
        print(f"  mutations on targets: {len(muts)}", file=sys.stderr)
    cnas: list[dict] = []
    if has_cna:
        cnas = fetch_cnas(profiles["COPY_NUMBER_ALTERATION"], study, list(entrez.values()))
        cnas_nonzero = [c for c in cnas if c.get("alteration") in CNA_MAP]
        print(f"  CNA events on targets: {len(cnas)} total ({len(cnas_nonzero)} non-diploid)", file=sys.stderr)
    svs: list[dict] = []
    if has_sv:
        svs = fetch_svs(profiles["STRUCTURAL_VARIANT"])
        print(f"  structural variants (study-wide): {len(svs)}", file=sys.stderr)

    # Bucket by sample
    by_mut: dict[str, list] = {}
    for m in muts:
        by_mut.setdefault(m["sampleId"], []).append(m)
    by_cna: dict[str, list] = {}
    for c in cnas:
        by_cna.setdefault(c["sampleId"], []).append(c)
    by_sv: dict[str, list] = {}
    for sv in svs:
        by_sv.setdefault(sv["sampleId"], []).append(sv)

    manifest_rows: list[dict] = []
    sample_dir = out_dir / study
    sample_dir.mkdir(parents=True, exist_ok=True)

    for s in samples:
        sid = s["sampleId"]
        m_list = by_mut.get(sid, [])
        c_list = by_cna.get(sid, [])
        sv_list = by_sv.get(sid, [])
        # Skip samples with no events on any target
        n_any_cna = sum(1 for c in c_list if c.get("alteration") in CNA_MAP)
        n_any_sv = sum(1 for sv in sv_list if sv.get("site1HugoSymbol") and sv.get("site2HugoSymbol") and sv["site1HugoSymbol"] != sv["site2HugoSymbol"])
        if not (m_list or n_any_cna or n_any_sv):
            continue

        clin = clinical.get(sid, {})
        subtype = classify_subtype(clin, sv_list)

        vcf_path = sample_dir / f"{sid}.vcf"
        cna_path = sample_dir / f"{sid}.cna.tsv"
        fusion_path = sample_dir / f"{sid}.fusion.tsv"
        n_mut = write_vcf(sid, m_list, vcf_path)
        n_cna = write_cna_tsv(sid, c_list, cna_path)
        n_sv = write_fusion_tsv(sid, sv_list, fusion_path)

        manifest_rows.append({
            "sample_id": sid,
            "study": study,
            "subtype": subtype,
            "pax_fusion": clin.get("PAX_FUSION", ""),
            "histology": clin.get("HISTOLOGICAL_SUBTYPE", ""),
            "risk_group": clin.get("RISK_GROUP", ""),
            "tumor_location": clin.get("PRIMARY_TUMOR_LOCATION", ""),
            "n_muts": n_mut,
            "n_cnas": n_cna,
            "n_fusions": n_sv,
            "vcf_path": str(vcf_path.relative_to(REPO_ROOT)),
            "cna_path": str(cna_path.relative_to(REPO_ROOT)),
            "fusion_path": str(fusion_path.relative_to(REPO_ROOT)),
        })
    print(f"  emitted {len(manifest_rows)} samples with at least one event", file=sys.stderr)
    return manifest_rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--studies", default="rms_nih_2014,rms_msk_2023",
                    help="Comma-separated cBioPortal study IDs.")
    ap.add_argument("--out-dir", type=Path,
                    default=REPO_ROOT / "data" / "target_rt")
    ap.add_argument("--manifest", type=Path,
                    default=REPO_ROOT / "data" / "target_rt" / "manifest.tsv")
    ap.add_argument("--targets-kb", type=Path,
                    default=REPO_ROOT / "assets" / "targets_kb.tsv")
    ap.add_argument("--gene-effect", type=Path, default=Path("/tmp/depmap/gene_effect.csv"))
    args = ap.parse_args()

    entrez = load_target_entrez(args.targets_kb, args.gene_effect)
    print(f"loaded Entrez IDs for {len(entrez)} target genes", file=sys.stderr)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict] = []
    for study in args.studies.split(","):
        all_rows.extend(fetch_one_study(study.strip(), entrez, args.out_dir))

    # Disambiguate sample IDs across studies if any collide
    seen: dict[str, int] = {}
    for r in all_rows:
        sid = r["sample_id"]
        if sid in seen:
            seen[sid] += 1
            r["sample_id"] = f"{sid}_{r['study']}"
        else:
            seen[sid] = 1

    cols = list(all_rows[0].keys()) if all_rows else []
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    with args.manifest.open("w", newline="") as fh:
        fh.write("## TARGET-RT + MSK-IMPACT cohort manifest (auto-generated by bin/fetch_target_rt.py)\n")
        fh.write(f"## Studies: {args.studies}\n")
        fh.write(f"## Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}\n")
        fh.write(f"## Samples written: {len(all_rows)}\n")
        if cols:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(all_rows)
    print(f"\nTOTAL: {len(all_rows)} samples written to {args.manifest}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
