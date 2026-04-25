# RMS-ISP Bundled Reference Data (v0.1)

These three TSVs are the **curated reference data** that lets the pipeline run end-to-end on a laptop in seconds without external API calls or large downloads. They are scoped to be just enough to demonstrate the pipeline's plumbing on the toy patient and the four pilot case studies.

| File | Purpose | Provenance | Replace before |
|------|---------|------------|----------------|
| `targets_kb.tsv` | Per-gene UniProt + known hotspot AA changes + LoF-as-driver flag, for the 21 RMS-relevant genes in `conf/targets.yaml` | Hand-curated by Mark Zobeck from v2 COG proposal target set, Shern 2014/2021, Crose 2012, Olanich 2015, standard RMS reviews | First scientific-claim run |
| `depmap_rms_summary.tsv` | Per-gene mean Chronos score across ~12 RMS cell lines, FP/FN-stratified, and a 0..1 dependency score | **Hand-curated PLACEHOLDER** consistent with published patterns (Dharia 2021 pediatric DepMap, Olanich CDK4, Marjon PRMT5/MTAP, Gryder BRD4) | Phase 3 real-data swap (DepMap 24Q2+ CRISPR_(DepMapPublic)\_v.csv filtered to RMS lineage) |
| `drug_target_map.tsv` | Curated drug-target pairs with mechanism, max clinical phase, and pediatric-evidence flag | Hand-curated from MTP, FDA pediatric labeling, COG/SARC trial registries (spot-checked, not exhaustive) | Phase 4 swap to live DGIdb v5 + OpenTargets + ClinicalTrials.gov |

## Why these are bundled rather than fetched

The pilot plan (§4 Aim 1, §11 Risk 1) is explicit that we should ship a runnable end-to-end pipeline before chasing data fidelity. These three files let the toy-patient integration test pass on any laptop with no network access, no API keys, and no cache downloads. They also document, in TSV form, exactly what shape the real upstream resources need to be reduced to before they enter the scoring engine; so the swap is a column-conformant file replacement, not a code change.

## What MUST be true before a real-data swap

1. The replacement file uses the same column names and types listed in the file's `## Schema:` header line.
2. The replacement file ships a `## PROVENANCE:` line identifying the source release (e.g., `DepMap 24Q2 from depmap.org/portal/download/all/, downloaded 2026-05-12, sha256 abc...`).
3. The replacement file is added to `/data/` (gitignored) and the path passed via `--depmap_summary` etc., never committed to the repo.

## Conditional drug applicability

Phase 4 currently honours one conditional rule encoded directly in the script: `KRAS_G12C_inhibitor` (adagrasib, sotorasib) is suppressed unless the matched variant is `G12C`. Future conditional rules (e.g., MDM2 inhibitors require TP53-wildtype tumor) belong in this same file as a future `applies_when` column rather than in the script.
