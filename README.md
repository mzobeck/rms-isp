# RMS-ISP: Rhabdomyosarcoma In Silico Pipeline

A modular, containerized Nextflow pipeline that takes molecular profiles from rhabdomyosarcoma tumors and outputs ranked, confidence-scored therapeutic hypotheses.

## Status

Pre-alpha. Pilot in progress. Not for clinical use.

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
nextflow run main.nf -profile laptop --input tests/data/toy_patient.vcf
```

This runs the toy RMS test patient end-to-end and produces a sample report. Used as the integration test in CI.

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

Part of the Zobeck Lab "in silico translational pipeline" program at Baylor College of Medicine / Texas Children's Hospital. See `../../proj_management/README.md` for the full project overview, aims, timeline, and `../../proj_management/PIPELINE_STANDUP_PLAN.md` for the step-by-step build plan.

## PI

Mark Zobeck, MD, MPH (BCM / Texas Children's Cancer and Hematology Center).
