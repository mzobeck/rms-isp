# MCI Transition Plan

Placeholder. To be completed once the public-data pilot reaches Aim 4 packaging.

## Intended contents

The exact configuration delta between the public-data pipeline and the MCI run. Items to document:

- Input file paths and dbGaP-authenticated access pattern.
- Reference genome build (both pilot and MCI use GRCh38; should be a no-op).
- Sample metadata schema: CCDI data model adapter (exercised against public CCDI submissions during pilot).
- QC thresholds: tumor purity, coverage depth (re-tuned on first 10 to 20 MCI samples).
- HIPAA-compliant compute environment (BCM RICC or institutional cloud).
- IRB protocol covering secondary analysis of de-identified dbGaP data.
- Regression test: every MCI run is compared to the public-data run on overlap samples; systematic differences flag pipeline or data problems.
