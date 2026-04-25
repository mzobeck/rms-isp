# Toy Test Patient

`toy_patient.vcf` is a hand-curated synthetic FN-RMS-leaning tumor designed to exercise the pipeline across multiple driver classes plus passenger negatives. It is not a real tumor and contains no PHI.

`toy_patient_expected.json` records, per variant, what the pipeline should produce at each phase. Used as the integration test harness.

## Contents

10 variants:

- 7 expected drivers spanning RAS-MAPK (NRAS Q61K, KRAS G12V, NF1 LoF), cell cycle / RB pathway (CDKN2A LoF), TP53 (R175H), MYOD1 (L122R neomorphic), FGFR4 (V550L kinase domain).
- 3 expected passengers: intronic, synonymous, high-frequency polymorphism.

## Coordinate status

Coordinates labeled `APPROX_verify_in_phase1` need verification against MANE/Ensembl 111 during Phase 1 implementation. The canonical hot-spots (TP53 R175H at chr17:7675088, NRAS Q61K at chr1:114713908, KRAS G12V at chr12:25245350) should be correct as written; if VEP fails to resolve them, that is a Phase 1 bug to investigate.

## How to use

```bash
nextflow run main.nf -profile laptop --input tests/data/toy_patient.vcf
```

The CI smoke test runs this on every commit. As phases come online, this same VCF is the integration-test input that exercises the full pipeline end-to-end.

## Disclosure

Synthetic. Variants are realistic in the sense that each gene-variant pair appears in the published RMS literature, but the combination of all of them in one tumor is artificial. Do not use for clinical or epidemiological inference.
