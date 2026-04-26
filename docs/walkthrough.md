# RMS-ISP Walkthrough

A guided tour of the Rhabdomyosarcoma In Silico Pipeline at v0.10.0-pilot. Read top-to-bottom for orientation; jump by section to dig in. Every claim links back to a file in the repo so you can ground out at code or data.

The pipeline takes a tumor's molecular profile (SNVs, copy-number alterations, fusions) and returns a ranked, confidence-scored list of candidate drugs with the full evidence trail behind each call. Five Nextflow phases run in Docker by default; output is one TSV per phase plus a clinician-readable markdown report. The bundled toy fixtures run end-to-end on a laptop in seconds, and the same pipeline runs on a 36-sample real-tumor cohort assembled from Shern 2014 TARGET-RT and MSK-IMPACT 2023.

## Table of contents

1. Why this pipeline exists
2. The five-phase architecture
3. The confidence formula
4. Running the pipeline
5. The case-study scorecard
6. Reference data: what is real, what is cached
7. Repository map
8. Roadmap

---

## 1. Why this pipeline exists

Rhabdomyosarcoma (RMS) is the most common pediatric soft-tissue sarcoma. The two molecular subtypes (fusion-positive PAX3/7::FOXO1 and fusion-negative) have distinct biology, but high-risk patients in both arms still die at unacceptable rates after standard VAC chemotherapy. Targeted therapy could move the needle, but only if a tumor's specific drivers are identified and matched to a drug with both mechanistic rationale and a plausible pediatric path.

RMS-ISP is the in-silico half of that workflow. Given a sequenced tumor, it returns a transparent, auditable ranking of drug hypotheses. Every score in the report links back to (a) the variant or event that triggered it, (b) the evidence database that supported it, and (c) the formula that combined them. The scientific gate is the four pilot case studies, codified as automated tests so the pipeline cannot silently regress.

The full design rationale lives in `plans/pilot_ccdi_mci/rms_translational_pilot_project.md` (gitignored, lives in the main repo dir).

## 2. The five-phase architecture

Each phase is a Nextflow process that wraps a single Python script in `bin/`. Inputs and outputs are TSVs; the contract between phases is "I read your TSV, I add columns, I write a new TSV." This makes every step independently testable and trivially debuggable: stop after any phase, open the TSV in a spreadsheet, see exactly what the next phase will see.

### Phase 1: variant + CNA + fusion annotation
`bin/phase1_annotate.py`, module `modules/phase1_variants/`

Input: a VCF, an optional CNA TSV, an optional fusion TSV. The annotator joins each event against the curated 21-gene target knowledge base in `assets/targets_kb.tsv` and emits one row per event with a `variant_call` of DRIVER, VUS, or PASSENGER and a numeric `variant_score` on 0 to 1. Hotspot residues, loss-of-function-bearing TSGs, oncogenic amplifications, and homozygous deletions of TSGs all classify as drivers; everything else falls to VUS or PASSENGER. v0.12 routes SNV annotation through a pluggable backend chain (`bin/annotators/`): toy fixtures use the `curated` backend (reads GENE / CONSEQUENCE from VCF INFO; scorecard stays fast and offline), real-cohort runs use `vep_rest` (Ensembl VEP REST API with an on-disk cache under `data/vep_cache/`; no 50 GB local install). v0.13 adds `oncokb` as a second-pass augmenting backend: when `ONCOKB_TOKEN` is set, the cohort runner upgrades to `vep_rest,oncokb`, and OncoKB's oncogenic call lifts non-hotspot missense variants to DRIVER. The pattern is documented in ADR 0003.

### Phase 2: structural reference attachment
`bin/phase2_structure.py`, module `modules/phase2_structure/`

Input: phase 1 output. For each event on a known target, attaches an AlphaFold structure ID and a `structural_score` reflecting how clean the structural rationale is for that event class (a hotspot residue in a kinase domain scores higher than a non-coding alteration on the same gene). v0.10 attaches references only; mutant-structure prediction with Boltz-1 / Chai-1 is gated to a v0.11 GPU profile run on a shortlist (see Roadmap).

### Phase 3: dependency integration
`bin/phase3_dependency.py`, module `modules/phase3_dependency/`

Input: phase 2 output, the DepMap RMS summary, and the OpenPedCan expression summary. Adds a subtype-aware `dependency_score` from real DepMap 26Q1 Chronos values across 14 RMS cell lines, plus an `expression_score` that captures RMS-vs-other-pediatric-tumor specificity from real OpenPedCan v15 RNA-seq (16 RMS samples, 4108 non-RMS pediatric samples). An oncogene-addiction floor lifts known driver dependencies (FGFR4 V550L, CDK4 amp, etc.) above what the cell-line average alone would suggest.

### Phase 4: drug matching
`bin/phase4_drugs.py`, module `modules/phase4_drugs/`

Input: phase 3 output, the curated drug-target map (`assets/drug_target_map.tsv`), the live DGIdb cache (`assets/dgidb_drugs.tsv`, ~195 rows over 21 genes after v0.16 noise filtering), and the live ClinicalTrials.gov cache (`assets/ctgov_rms_drugs.tsv`). Joins each annotated event to all drugs whose targets match. Each drug-event pair carries a `drug_evidence_score` that combines phase-of-development weight with pediatric-evidence weight (yes_approved > pediatric_trial > adult_only > preclinical). Output is long-format: one row per variant-drug pair, so a single hotspot mutation can fan out into many candidate drugs. v0.16 added a noise filter in `bin/fetch_dgidb.py` that drops interactions whose only source is a clinical-trial-listing database (e.g. ClearityFoundationClinicalTrial), since those typically reflect drug + biomarker co-occurrence in trial enrollment criteria rather than mechanistic targeting.

### Phase 5: confidence scoring and report
`bin/phase5_score.py`, module `modules/phase5_scoring/`

Combines the four upstream component scores into a single `confidence` value per variant-drug pair using the fixed weights below, then writes both a long TSV and a human-readable markdown report at `results/phase5/<sample>.report.md`. The report header records pipeline version, input file SHA, per-component weights, and the formula itself so a reviewer can audit any score.

## 3. The confidence formula

A row's confidence is a fixed-weight linear combination of five component scores, each on 0 to 1:

| Component  | Source                              | Weight |
|------------|-------------------------------------|--------|
| variant    | phase 1 (`variant_score`)           | 0.25   |
| structural | phase 2 (`structural_score`)        | 0.15   |
| dependency | phase 3 (`dependency_score`)        | 0.25   |
| expression | phase 3 (`expression_score`)        | 0.15   |
| drug       | phase 4 (`drug_evidence_score`)     | 0.20   |

`confidence = sum over k of (weight_k * component_k)`. Weights sum to 1.0; the formula is defined in pilot plan §4 Aim 2 Phase 5 and copied verbatim into `bin/phase5_score.py`.

As of v0.10 every component draws from live data: DepMap for dependency, OpenPedCan for expression, DGIdb for drug interactions, ClinicalTrials.gov for pediatric-evidence weighting. v0.12 added live VEP annotation for the variant component on real cohorts; v0.13 added an OncoKB augmentation pass that lifts non-hotspot missense variants to DRIVER when OncoKB calls them oncogenic. Toys keep curated INFO for scorecard speed. The structural component is still curated rather than computed at scale; Boltz/Chai mutant-structure prediction is the GPU-gated future priority.

The rigid scoring rule is part of the design. If a case study breaks, the fix is to repair the broken phase or its reference data, never to retune the weights. That discipline is recorded as Risk 6 in the pilot plan; v0.5 added `toy_fgfr4_only.vcf` and v0.9 added `toy_ras_mek.vcf` for exactly this reason (when an assertion failed, the response was a new fixture that pinned the broken behavior, not a weight change).

## 4. Running the pipeline

The default invocation runs the bundled toy patient through all five phases inside Docker:

```bash
containers/build.sh                  # one-time: build rms-isp/base:0.10.0
nextflow run main.nf -profile laptop
```

Without a Docker daemon, fall back to host Python:

```bash
nextflow run main.nf -profile laptop_nodocker
```

A different fixture and subtype:

```bash
nextflow run main.nf -profile laptop \
    --input tests/data/toy_fp_cdk4amp.vcf \
    --cna tests/data/toy_fp_cdk4amp.cna.tsv \
    --fusion tests/data/toy_fp_cdk4amp.fusion.tsv \
    --sample_id TOY_FP_CDK4 \
    --subtype FP
```

Run all three toy patients plus a cross-sample summary:

```bash
bin/run_all_toys.sh   # writes results/multisample_summary.md
```

Run the 36-sample real-tumor cohort (Shern 2014 TARGET-RT + MSK-IMPACT 2023):

```bash
python3 bin/run_target_rt.py   # writes results/target_rt_cohort_summary.tsv
```

Outputs land in `results/<phaseN>/<sample_id>.<phaseN>.tsv`, with the final markdown report at `results/phase5/<sample_id>.report.md`. Nextflow trace, timeline, and DAG go to `results/pipeline_info/`.

Param defaults are in `nextflow.config`; profile configs (laptop, laptop_nodocker, slurm, aws, gpu) are in `conf/`.

## 5. The case-study scorecard

Pilot acceptance is encoded as ten assertions across five cases in `tests/cases.toml` and run by `bin/check_case_studies.py`. CI runs the scorecard on every push and blocks merge if any assertion regresses.

| #       | Case                                  | Fixture        | Asserts |
|---------|---------------------------------------|----------------|---------|
| 1       | MTAP / CDKN2A 9p21 co-deletion        | TOY_MTAP_NULL  | PRMT5 inhibitors top 7, CDK4/6 inhibitors top 3 |
| 2       | CDK4 amplification in FP-RMS          | TOY_FP_CDK4    | CDK4/6 inhibitors top 3, CDK4 in top 3, BET inhibitors top 6 |
| 3       | FGFR4 V550L hotspot                   | TOY_FGFR4      | FGFR inhibitors top 3, FGFR4 in top 3 |
| 4       | NRAS Q61K + MEK                       | TOY_RAS_MEK    | MEK inhibitors top 2, NRAS in top 2 |
| general | passenger sanity check                | TOY_TUMOR      | every PASSENGER row below 0.10 confidence |

Three assertion kinds:

- `mechanism_in_top_n`: at least one drug in the top N has a mechanism string containing the substring (e.g. "PRMT5_inhibitor").
- `gene_event_in_top_n`: the named gene appears in the top N rows.
- `passenger_below`: every PASSENGER-classified row scores below the threshold.

A PASS means the pipeline reproduces the scientific call we know to be correct from the literature. It does not mean the score itself is calibrated to a clinical decision threshold; calibration is a v1.0 concern.

```bash
python3 bin/check_case_studies.py    # runs all 5 cases, all 10 asserts
open results/scorecard.md
```

## 6. Reference data: what is real, what is cached

Everything in `assets/` is one of three kinds:

- **Hand-curated, versioned with the repo.** `targets_kb.tsv` (the 21-gene knowledge base), `drug_target_map.tsv` (the curated drug catalog), `rms_cell_lines.tsv` (the DepMap cell-line registry).
- **Live cache, refreshable, versioned with the repo.** `dgidb_drugs.tsv`, `ctgov_rms_drugs.tsv`, `depmap_rms_summary.tsv`, `openpedcan_expression_summary.tsv`. Each was generated by a `bin/fetch_*.py` script and carries a header comment with source URL and generation timestamp. CI does not refresh these; rerunning the fetcher and committing the new TSV is a deliberate choice.
- **Empty placeholders.** `empty.cna.tsv` and `empty.fusion.tsv` so SNV-only patients still satisfy the channel signatures.

Refreshing a live cache:

```bash
python3 bin/fetch_dgidb.py             # 21 genes  -> assets/dgidb_drugs.tsv
python3 bin/fetch_clinicaltrials.py    # ClinicalTrials.gov pediatric trials
python3 bin/fetch_depmap.py            # DepMap 26Q1 Chronos, 14 RMS lines
/tmp/openpedcan_venv/bin/python bin/fetch_openpedcan_expression.py
python3 bin/fetch_target_rt.py         # 36-sample real-tumor cohort
```

The OpenPedCan fetcher needs its own venv because it depends on `pyreadr` to read the upstream `.rds` file. The pipeline itself is stdlib-only; the venv is fetcher-only.

## 7. Repository map

```
rms-isp/
├── main.nf                Nextflow workflow: 5 phases wired together
├── nextflow.config        defaults + manifest (version, params)
├── conf/                  per-profile configs (laptop, slurm, aws, gpu)
├── modules/
│   ├── phase1_variants/   Nextflow process wrappers, one per phase
│   ├── phase2_structure/
│   ├── phase3_dependency/
│   ├── phase4_drugs/
│   └── phase5_scoring/
├── bin/                   Python that does the science (one script per phase + fetchers)
├── assets/                reference data the science consumes
├── containers/            base.Dockerfile + 5 per-phase Dockerfiles + build.sh
├── tests/
│   ├── data/              toy patients (VCF, CNA, fusion) + golden expected outputs
│   └── cases.toml         scorecard assertions (source of truth)
└── docs/                  this walkthrough + architecture + data sources + ADRs
```

Three things to know when navigating:

- `bin/` is where the actual computation lives; `modules/` is thin Nextflow process glue.
- TSVs in `assets/` with a comment-block header are auto-generated; do not hand-edit.
- The `plans/` folder lives in the main repo dir and is gitignored. Only the build plan and the data overview are canonical; the rest is scratch.

## 8. Roadmap

In rough priority order:

1. **Real VEP + OncoKB for phase 1** (shipped v0.12 + v0.13 via REST API). `bin/annotators/{vep_rest,oncokb}.py` form a chained annotator: VEP first, OncoKB second. Per-variant responses cache to `data/vep_cache/` and `data/oncokb_cache/`. OncoKB requires an academic API token (`ONCOKB_TOKEN` env var); without it, the chain falls back to VEP-only. Local-CLI VEP (50 GB cache, REVEL / AlphaMissense plugins) remains a future option behind a separate backend if rate limits or plugin needs ever bite.
2. **Cohort-level visualizations** (shipped v0.11). Static SVG charts embedded in `results/target_rt_cohort_summary.md`: mechanism distribution, per-target druggability (gene by subtype, scales to any N), and a per-sample heatmap rendered for cohorts up to 100 samples. Implementation in `bin/cohort_visualize.py`.
3. **Tier classification in phase 5** (shipped v0.14 per-row + v0.15 portfolio). Per-row tiers in `bin/phase5_score.py` (Tier 1 = FDA-approved drug + DRIVER variant; Tier 2 = phase 2/3 + DRIVER; Tier 3 = phase 1/preclinical + DRIVER; "" = non-DRIVER). v0.15 adds cohort-aggregated promotion: `bin/cohort_portfolio.py` computes per-(gene, subtype) prevalence, picks the best-available drug per gene, and applies the v2 plan §5 prevalence gates (Tier 1 = approved + >=5% of any subtype; Tier 2 = late-phase + >=3%; Tier 3 = qualifying DRIVER hits that did not promote). The output is `results/target_rt_STS_committee_portfolio.md`, the deliverable referenced in v2 plan §5.
4. **Live OpenTargets in phase 4.** Was deferred at v0.6 because the v4 GraphQL schema does not expose drugs directly on `Target`. Adds disease-association evidence to drug ranking.
5. **Boltz-1 / Chai-1 mutant prediction in phase 2.** GPU-bound, gates on the `gpu` profile. Run on a shortlist of top-N target-drug pairs per tumor rather than every event.

Eventually: the COG Molecular Characterization Initiative handoff (pilot plan §12), which is mostly a CCDI-to-internal-schema adapter plus dbGaP-authenticated compute and IRB.

## Where to look for more

- `README.md` for the elevator pitch and quickstart.
- `docs/architecture.md`, `docs/data_sources.md`, `docs/MCI_TRANSITION.md` for component-level detail.
- `docs/adr/` for the architectural decisions and their rationale.
- `tests/cases.toml` for the scorecard contract.
- `bin/phase5_score.py` for the confidence formula in code.
- `plans/pilot_ccdi_mci/rms_translational_pilot_project.md` (gitignored) for the canonical build plan.
