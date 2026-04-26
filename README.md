# RMS-ISP: Rhabdomyosarcoma In Silico Pipeline

A modular, containerized Nextflow pipeline that takes molecular profiles from rhabdomyosarcoma tumors and outputs ranked, confidence-scored therapeutic hypotheses.

## Status

**v0.8.0-pilot**: Real-tumor cohort now **36 samples** (Shern 2014 TARGET-RT + MSK-IMPACT 2023) running on **SNV + CNA + fusion** data from cBioPortal. FP-RMS samples surface via PAX3-FOXO1 / PAX7-FOXO1 fusion calls that auto-route to BET inhibitors (14 samples). CDK4 / MDM2 amplifications surface via CNA calls and route to CDK4/6 inhibitors (3 samples). The cohort-level mechanism distribution recapitulates the textbook RMS subtype-to-therapy logic without any retuning: BET inhibitors for FP, MEK for FN RAS-MAPK, CDK4/6 for amplifications, FGFR for FGFR4 hotspots, PI3K for PIK3CA, WEE1 for TP53. All five toy case studies still PASS (9/9 assertions). Phase 3 uses real DepMap 26Q1 + oncogene-addiction floor. Phase 4 unions curated drug map with live DGIdb + ClinicalTrials.gov upgrades. Not for clinical use.

## Mission

Convert standardized pediatric cancer molecular data (TARGET-RT, OpenPedCan, DepMap, and eventually the COG Molecular Characterization Initiative) into ranked, evidence-scored therapeutic hypotheses ready for experimental validation.

The pipeline is designed to be:

- Reproducible: every run records pipeline version, container digests, parameters, and input data version.
- Modular: each phase is an independently testable Nextflow process.
- Portable: runs on laptop, institutional HPC, or cloud (AWS Batch, GCP Batch).
- Generalizable: phases 1 to 5 are disease-agnostic; RMS is the first application.

## Pipeline phases

1. Variant discovery and annotation (VEP, ANNOVAR, CADD, REVEL, AlphaMissense, OncoKB, ClinVar, COSMIC).
2. Structural modeling and druggability (AlphaFold DB, FoldX, Boltz-1/Chai-1 for shortlist, AutoDock Vina).
3. Expression and dependency integration (DepMap Chronos, OpenPedCan RNA-seq, STRING, clusterProfiler).
4. Drug matching (DGIdb, OpenTargets, Molecular Targets Platform, ChEMBL, DrugBank, ClinicalTrials.gov, CMap/LINCS L1000).
5. Confidence scoring and patient report generation.

Each phase ships with positive and negative controls and a defined validation gate.

## Quickstart

```bash
# default: TOY_TUMOR (FN-RMS, SNV-only, exercises case studies 3 + 4)
nextflow run main.nf -profile laptop

# CNA + fusion example: TOY_FP_CDK4 (FP-RMS, exercises case study 2)
nextflow run main.nf -profile laptop \
    --input tests/data/toy_fp_cdk4amp.vcf \
    --cna tests/data/toy_fp_cdk4amp.cna.tsv \
    --fusion tests/data/toy_fp_cdk4amp.fusion.tsv \
    --sample_id TOY_FP_CDK4 \
    --subtype FP

# Run all three bundled toy patients + cross-sample summary:
bin/run_all_toys.sh   # writes results/multisample_summary.md
```

Outputs land in `results/`:

```
results/
├── phase1/TOY_TUMOR.phase1.tsv         annotated variants (DRIVER / VUS / PASSENGER calls)
├── phase2/TOY_TUMOR.phase2.tsv         + AlphaFold structural reference + structural score
├── phase3/TOY_TUMOR.phase3.tsv         + DepMap dependency score (subtype-aware)
├── phase4/TOY_TUMOR.phase4.tsv         + matched drugs (long format, one row per variant-drug)
├── phase5/
│   ├── TOY_TUMOR.phase5.tsv            + confidence score and per-component scores
│   └── TOY_TUMOR.report.md             ranked therapeutic-hypothesis report (human-readable)
└── pipeline_info/                      Nextflow trace, timeline, report, dag
```

### Case-study scorecard

Pilot acceptance criteria (per `plans/pilot_ccdi_mci/rms_translational_pilot_project.md` §3) are encoded in `tests/cases.toml` and asserted by `bin/check_case_studies.py`. CI runs the scorecard on every push and blocks merge if any assertion regresses.

```bash
python3 bin/check_case_studies.py
open results/scorecard.md
```

| Pilot # | Test fixture | Asserts |
|---|---|---|
| 1 (MTAP/PRMT5) | TOY_MTAP_NULL | PRMT5 inhibitors in top 7, CDK4/6 inhibitors in top 3 |
| 2 (CDK4 amp) | TOY_FP_CDK4 | CDK4/6 inhibitors in top 3, CDK4 in top 3, BET inhibitors in top 6 |
| 3 (FGFR4) | TOY_TUMOR | FGFR inhibitors in top 5 |
| 4 (RAS/MEK) | TOY_TUMOR | MEK inhibitors at #1-2 |
| general | TOY_TUMOR | All 3 passengers below 0.10 confidence |

## Repository layout

```
rms-isp/
├── main.nf                         Top-level workflow
├── nextflow.config                 Default config
├── conf/                           Profile configs (laptop, slurm, aws, gpu)
├── modules/                        One subdirectory per pipeline phase
├── containers/                     Dockerfile per module
├── tests/
│   ├── data/                       Toy patient + golden outputs
│   └── integration/                End-to-end tests
├── docs/                           Architecture, data sources, MCI transition
└── .github/workflows/              CI
```

## Documentation

- `docs/architecture.md`: pipeline architecture (placeholder)
- `docs/data_sources.md`: public data sources, manifests, licenses (placeholder)
- `docs/MCI_TRANSITION.md`: configuration delta for MCI data (placeholder)
- `docs/adr/`: architecture decision records

## License

Code: MIT.
Documentation and reference outputs: CC-BY-4.0.
Data sources retain their original licenses; see `LICENSES.md`.

## Project context

Part of the Zobeck Lab "in silico translational pipeline" program at Baylor College of Medicine / Texas Children's Hospital. The canonical build plan is `plans/pilot_ccdi_mci/rms_translational_pilot_project.md`; the surrounding data landscape is in `plans/pilot_ccdi_mci/ccdi_mci_data_overview.md`. The `plans/` folder is gitignored.

## PI

Mark Zobeck, MD, MPH (BCM / Texas Children's Cancer and Hematology Center).
