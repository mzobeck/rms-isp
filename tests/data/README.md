# Toy Test Patients

Three hand-curated synthetic tumors that, taken together, exercise the four pilot case studies from `plans/pilot_ccdi_mci/rms_translational_pilot_project.md` §3. None are real tumors; none contain PHI.

| Sample | Profile | Pilot case study | Files |
|---|---|---|---|
| `TOY_TUMOR` | FN-RMS, 7 driver SNVs + 3 passengers | 3 (FGFR4) and 4 (RAS/MEK) | `toy_patient.vcf` |
| `TOY_FP_CDK4` | FP-RMS, PAX3-FOXO1 fusion + CDK4/MDM2 amplification + TP53 R175H | 2 (CDK4 amp → CDK4/6i) | `toy_fp_cdk4amp.{vcf,cna.tsv,fusion.tsv}` |
| `TOY_MTAP_NULL` | FN-RMS, CDKN2A + MTAP 9p21 homozygous co-deletion | 1 (MTAP/CDKN2A → PRMT5 + CDK4/6i) | `toy_mtap_null.{vcf,cna.tsv,fusion.tsv}` |

`toy_patient_expected.json` records, per variant, what the pipeline should produce at each phase for the original SNV-only fixture. Used as the integration test harness.

## Contents

10 variants:

- 7 expected drivers spanning RAS-MAPK (NRAS Q61K, KRAS G12V, NF1 LoF), cell cycle / RB pathway (CDKN2A LoF), TP53 (R175H), MYOD1 (L122R neomorphic), FGFR4 (V550L kinase domain).
- 3 expected passengers: intronic, synonymous, high-frequency polymorphism.

## Coordinate status

Coordinates labeled `APPROX_verify_in_phase1` need verification against MANE/Ensembl 111 during Phase 1 implementation. The canonical hot-spots (TP53 R175H at chr17:7675088, NRAS Q61K at chr1:114713908, KRAS G12V at chr12:25245350) should be correct as written; if VEP fails to resolve them, that is a Phase 1 bug to investigate.

## How to use

Run a single sample:

```bash
nextflow run main.nf -profile laptop \
    --input tests/data/toy_fp_cdk4amp.vcf \
    --cna tests/data/toy_fp_cdk4amp.cna.tsv \
    --fusion tests/data/toy_fp_cdk4amp.fusion.tsv \
    --sample_id TOY_FP_CDK4 \
    --subtype FP
```

Run all three toy patients with a cross-sample summary:

```bash
bin/run_all_toys.sh
open results/multisample_summary.md
```

The CI smoke test runs `TOY_TUMOR` on every commit. The other two patients exist to keep pilot case studies 1 and 2 covered as Phase 1 and Phase 4 evolve.

## Disclosure

Synthetic. Variants are realistic in the sense that each gene-variant pair appears in the published RMS literature, but the combination of all of them in one tumor is artificial. Do not use for clinical or epidemiological inference.
