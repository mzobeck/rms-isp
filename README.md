# RMS-ISP: Rhabdomyosarcoma In Silico Pipeline

A modular, containerized Nextflow pipeline that takes molecular profiles from rhabdomyosarcoma tumors and outputs ranked, confidence-scored therapeutic hypotheses.

## Status

**v0.2.0-pilot** — end-to-end with SNV + CNA + fusion inputs. Three toy patients exercising all four pilot case studies (FGFR4, RAS/MEK, CDK4 amp, MTAP/PRMT5) all land their expected drug classes near the top. Runs in seconds on a laptop, no external API calls or large downloads required. Reference data for Phase 3 (DepMap) and Phase 4 (drug-target map) is curated and bundled in `assets/`; both swap to live data sources without code changes (see `assets/README.md`). Not for clinical use.

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

### v0.2 eyeball-tests (one per case study)

| Case | Test fixture | Expected top result |
|---|---|---|
| 1 (MTAP/PRMT5) | TOY_MTAP_NULL | CDK4/6 inhibitors via CDKN2A loss in top 3 + PRMT5 inhibitors (MRTX1719 class) via MTAP loss in top 7 |
| 2 (CDK4 amp) | TOY_FP_CDK4 | palbociclib / ribociclib / abemaciclib at #1-3; PAX3-FOXO1 → BET inhibitors in top 6; MDM2 → MDM2 inhibitors close behind |
| 3 (FGFR4) | TOY_TUMOR | erdafitinib + futibatinib in top 5 |
| 4 (RAS/MEK) | TOY_TUMOR | trametinib + selumetinib at #1-2 |

All three passenger variants in TOY_TUMOR (TP53 intronic, NRAS/KRAS synonymous) should rank at the bottom of their respective reports.

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
