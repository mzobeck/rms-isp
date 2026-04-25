# RMS-ISP Walkthrough

A guided tour of the Rhabdomyosarcoma In Silico Pipeline (RMS-ISP). Designed for voice-chat use: a Claude can read this top-to-bottom for a high-level pass, then zoom into any section by anchor when you want detail. Every claim links back to a specific file or external reference so you can ground out at code, data, or literature.

The pipeline takes a tumor's molecular profile (variants + copy-number alterations + fusions) and returns a ranked, confidence-scored list of candidate drugs with the full evidence trail behind each call. v0.4 runs end-to-end on a laptop in seconds against three bundled toy patients that exercise all four pilot acceptance case studies.

## How to use this document

Three reading paths:

1. **High-level pass**: read [Section 1](#1-the-problem-and-the-design) and [Section 2](#2-the-architecture-in-one-page), then skim [Section 3](#3-tour-of-the-repository). About 10 minutes spoken.
2. **Phase-by-phase deep dive**: jump to [Section 5](#5-the-five-phases-in-depth). Each phase has the same shape: inputs, what it does, outputs, the algorithm or knowledge it depends on, and what's deferred. About 30 minutes spoken.
3. **Concept zoom**: when a biology or algorithm term shows up, jump to [Appendix A: Glossary](#appendix-a-glossary-of-biological-and-computational-concepts). Each glossary entry stands alone.

## Table of contents

1. [The problem and the design](#1-the-problem-and-the-design)
2. [The architecture in one page](#2-the-architecture-in-one-page)
3. [Tour of the repository](#3-tour-of-the-repository)
4. [The reference data, in detail](#4-the-reference-data-in-detail)
5. [The five phases, in depth](#5-the-five-phases-in-depth)
   - [5.1 Phase 1: variant + CNA + fusion annotation](#51-phase-1-variant--cna--fusion-annotation)
   - [5.2 Phase 2: structural reference attachment](#52-phase-2-structural-reference-attachment)
   - [5.3 Phase 3: dependency integration](#53-phase-3-dependency-integration)
   - [5.4 Phase 4: drug matching](#54-phase-4-drug-matching)
   - [5.5 Phase 5: confidence scoring and report](#55-phase-5-confidence-scoring-and-report)
6. [End-to-end: what happens when you press go](#6-end-to-end-what-happens-when-you-press-go)
7. [How we know it works: the case-study scorecard](#7-how-we-know-it-works-the-case-study-scorecard)
8. [What's real, what's placeholder, what's deferred](#8-whats-real-whats-placeholder-whats-deferred)
9. [Where the pipeline goes from here](#9-where-the-pipeline-goes-from-here)
10. [Appendix A: glossary of biological and computational concepts](#appendix-a-glossary-of-biological-and-computational-concepts)
11. [Appendix B: file-by-file reference](#appendix-b-file-by-file-reference)

---

## 1. The problem and the design

### 1.1 The clinical problem

Rhabdomyosarcoma (RMS) is the most common pediatric soft-tissue sarcoma. Despite roughly two decades of detailed molecular characterization (IRS-IV through ARST1431), overall survival for high-risk RMS has not improved since the 1990s. The bottleneck is not data: there are large consortium cohorts (INSTRuCT, around 6,972 patients), deep epigenetic maps (Gryder 2017-2020), single-cell atlases (Wei 2022, Patel 2022), tumor microenvironment work (DeMartino 2023), and resistant-state transcriptomics (Danielli 2024). The bottleneck is the absence of a disciplined pipeline that takes those data and returns actionable drug hypotheses at the level of an individual tumor.

The COG Molecular Characterization Initiative (MCI, dbGaP `phs002790`) is generating standardized whole-exome, fusion-panel, and methylation data on every consenting newly-diagnosed pediatric soft-tissue sarcoma in the United States. Within 2-3 years the MCI will have over 500 RMS cases, the largest standardized molecular dataset for the disease. The pipeline this repo builds is the engineering substrate that will turn those cases into ranked therapeutic-hypothesis reports for the COG Soft Tissue Sarcoma committee.

### 1.2 The design philosophy

Five principles, all inherited from `plans/pilot_ccdi_mci/rms_translational_pilot_project.md`:

1. **Build the pipeline on public data first**. Validate every phase against a literature ground truth before MCI data ever touches it. When MCI access arrives, swapping the input is a configuration change, not a rewrite. We are not debugging the pipeline at MCI handoff, we are interpreting results.
2. **Modular and independently testable**. Each phase is a separate Nextflow process with a defined input and output schema. Any phase can be re-implemented (or replaced with a real version of a placeholder) without touching downstream phases.
3. **Auditable scoring**. The confidence formula is a simple weighted sum with fixed weights from the pilot plan, and every component is exposed in the per-event detail of the report. A clinician or reviewer can see exactly why a drug ranked where it did.
4. **Reproducible by construction**. Every run records pipeline version, container digests (when containerized), parameters, and input data version. Toy fixtures are committed to the repo so the integration test runs cold on any machine.
5. **No retrofitting**. From the pilot plan §11 Risk 6 and §15: if a case study returns the wrong top drug, fix the broken phase, do not widen the confidence weights to rescue the expected answer. The case-study scorecard exists to catch this.

### 1.3 What the pipeline produces

For each input tumor: a ranked Markdown report listing the top target-drug pairs with the full evidence trail behind each call. The header includes pipeline version, input file SHA, run timestamp, and a v0.x disclaimer. The body has a top-10 table, a per-event detail block (one per variant / CNA / fusion) that exposes every component score and links to the AlphaFold reference structure, and a methodology section that prints the scoring weights in plain text.

A real example lives at `results/phase5/TOY_TUMOR.report.md` after running the default pipeline. See [Section 6](#6-end-to-end-what-happens-when-you-press-go) for how to produce it.

---

## 2. The architecture in one page

```
                  +--------------------+
                  |  inputs per tumor  |
                  |  VCF + CNA + Fusion|
                  +---------+----------+
                            |
                            v
+---------------------------+----------------------------+
| Phase 1: variant / CNA / fusion annotation             |
| in:  VCF, CNA TSV, fusion TSV, targets_kb.tsv          |
| out: per-event TSV with DRIVER / VUS / PASSENGER call  |
+---------------------------+----------------------------+
                            |
                            v
+---------------------------+----------------------------+
| Phase 2: structural reference                          |
| in:  Phase 1 TSV                                        |
| out: + AlphaFold URL + structural_score                 |
+---------------------------+----------------------------+
                            |
                            v
+---------------------------+----------------------------+
| Phase 3: dependency integration                        |
| in:  Phase 2 TSV, depmap_rms_summary.tsv               |
| out: + dependency_score (subtype-aware)                 |
+---------------------------+----------------------------+
                            |
                            v
+---------------------------+----------------------------+
| Phase 4: drug matching                                 |
| in:  Phase 3 TSV, drug_target_map.tsv, dgidb_drugs.tsv |
| out: long-format TSV (one row per event-drug pair)     |
+---------------------------+----------------------------+
                            |
                            v
+---------------------------+----------------------------+
| Phase 5: confidence scoring + report                   |
| in:  Phase 4 TSV, original VCF                         |
| out: scored TSV + ranked Markdown report               |
+---------------------------+----------------------------+
                            |
                            v
                  +---------+----------+
                  |  results/phase5/   |
                  |  <sample>.report.md|
                  +--------------------+
```

Each box is a Nextflow process backed by a Python script. The arrow between them is a Nextflow channel carrying a tuple of (sample_id, file path). The boxes can run on a laptop, a Slurm cluster, AWS Batch, or GCP Batch by switching the `-profile` flag.

### 2.1 Why Nextflow

[ADR 0002](adr/0002-nextflow-orchestrator.md). The shortlist was Nextflow vs Snakemake. Nextflow won on:

- Native first-class support for AWS Batch, GCP Batch, and Slurm without a rewrite. Same code, different `-profile`.
- nf-core community has production-grade module patterns directly applicable to genomic workflows.
- PedcBioPortal, OpenPedCan, and most CCDI-affiliated pipelines use Nextflow. Cross-pipeline reuse is cheaper inside that ecosystem.
- Seqera Platform (Tower) is a drop-in for provenance and monitoring later.

We commit to DSL2 modules. Process-level idioms follow nf-core conventions where applicable.

### 2.2 Why Python scripts behind every Nextflow process

Each `modules/phaseN_<name>/main.nf` is a thin Nextflow wrapper that invokes a Python script in `bin/`. The Python script does the actual work. Two reasons:

1. The Python script is independently testable from the command line, no Nextflow required. This makes development fast and CI cheap (the case-study scorecard at [bin/check_case_studies.py](../bin/check_case_studies.py) runs the Python directly).
2. The Nextflow layer adds orchestration, provenance, containerization, and resource management on top, without any of those concerns leaking into the science code.

Rule: Python contains the science. Nextflow contains the plumbing.

### 2.3 Why TSVs everywhere

Every phase reads and writes tab-separated values with a documented schema. Reasons:

- Human-readable. You can `cat` any intermediate output and see what's going on. The pilot plan §11 Risk 5 calls out that interpretability matters more than computational elegance for v0.x.
- Tooling-agnostic. Works with pandas, polars, R, awk, sqlite import, anything.
- Diff-able under git when small.
- Schema is enforced by the Python `csv.DictReader` + explicit column lists. Any schema change forces an update to every downstream phase, which is desirable.

When phase outputs grow past a few thousand rows we will move to Parquet for the on-disk format while keeping the column contract identical. That is a v1.0 change, not now.

---

## 3. Tour of the repository

Top-level layout, with the parts you will touch most marked:

```
rms-isp/
├── README.md                        ← project overview, quickstart, eyeball-tests
├── main.nf                          ← top-level Nextflow workflow
├── nextflow.config                  ← all param defaults; CLI overrides
├── conf/
│   ├── base.config                  ← shared resource defaults
│   ├── laptop.config                ← local CPU profile
│   ├── slurm.config                 ← institutional HPC
│   ├── aws.config                   ← AWS Batch
│   ├── gpu.config                   ← GPU nodes for future Boltz/Chai work
│   └── targets.yaml                 ← curated list of 21 RMS-relevant target genes
├── modules/                         ← one Nextflow module per phase
│   ├── phase1_variants/main.nf
│   ├── phase2_structure/main.nf
│   ├── phase3_dependency/main.nf
│   ├── phase4_drugs/main.nf
│   └── phase5_scoring/main.nf
├── bin/                             ← Python implementations
│   ├── phase1_annotate.py
│   ├── phase2_structure.py
│   ├── phase3_dependency.py
│   ├── phase4_drugs.py
│   ├── phase5_score.py
│   ├── fetch_dgidb.py               ← refresh DGIdb cache
│   ├── check_case_studies.py        ← scorecard runner (CI uses this)
│   └── run_all_toys.sh              ← runs all 3 toy patients + scorecard
├── assets/                          ← bundled reference data (in-repo, small)
│   ├── targets_kb.tsv               ← per-gene biology (hotspots, LoF flags, amp/del flags)
│   ├── depmap_rms_summary.tsv       ← Chronos summary, PLACEHOLDER
│   ├── drug_target_map.tsv          ← curated drug-target map
│   ├── dgidb_drugs.tsv              ← auto-generated DGIdb cache
│   ├── empty.cna.tsv                ← default empty CNA file
│   └── empty.fusion.tsv             ← default empty fusion file
├── tests/
│   ├── cases.toml                   ← acceptance criteria for the scorecard
│   └── data/
│       ├── toy_patient.vcf          ← TOY_TUMOR (FN-RMS, SNV-only, 7 drivers + 3 passengers)
│       ├── toy_patient_expected.json
│       ├── toy_fp_cdk4amp.{vcf,cna.tsv,fusion.tsv}    ← TOY_FP_CDK4 (case study 2)
│       └── toy_mtap_null.{vcf,cna.tsv,fusion.tsv}     ← TOY_MTAP_NULL (case study 1)
├── containers/                      ← Dockerfile per phase (placeholders)
├── docs/
│   ├── walkthrough.md               ← this file
│   ├── architecture.md              ← short note pointing at the pilot plan
│   ├── data_sources.md              ← placeholder; populated at Aim 1 close
│   ├── MCI_TRANSITION.md            ← placeholder; the eventual config delta
│   └── adr/                         ← architecture decision records
├── plans/                           ← gitignored; lives in the main repo dir only
│   └── pilot_ccdi_mci/
│       ├── rms_translational_pilot_project.md   ← canonical build plan
│       └── ccdi_mci_data_overview.md            ← data landscape reference
└── .github/workflows/ci.yml         ← scorecard + Nextflow + markdown lint
```

Three folders carry the load:

- `bin/`: the Python that does the science.
- `assets/`: the reference data the science consumes.
- `tests/`: the toy patients + the scorecard's machine-checkable acceptance criteria.

Everything else is orchestration, config, documentation, or future scaffolding.

### 3.1 The plans folder is gitignored on purpose

The pilot plan and data overview live in `plans/` in the main repo dir, NOT in the worktree, NOT in git. Reason: the source documents are .docx files (large binaries) plus markdown exports, and they are working drafts that change frequently. Treating them as living scientific documents outside the code repo is intentional. The README and `docs/architecture.md` reference the canonical paths.

---

## 4. The reference data, in detail

Three TSVs in `assets/` are what makes the pipeline interpretable. Each one documents in its header what it is, where it came from, and what must replace it before any scientific claim is made.

### 4.1 `targets_kb.tsv`: the per-gene knowledge base

[assets/targets_kb.tsv](../assets/targets_kb.tsv). 21 genes, hand-curated by Mark Zobeck from the v2 COG proposal target set. Schema (one row per gene):

| Column | What it does |
|---|---|
| `gene` | HUGO symbol |
| `uniprot` | Reference UniProt accession (drives Phase 2 AlphaFold lookup) |
| `role` | `oncogene` / `tsg` / `tf` / `enzyme` / `epigenetic` / `metabolic` / `fusion_partner` |
| `hotspots` | Comma-separated AA changes (HGVSp short, e.g., `Q61K,Q61L,G12V`); a missense variant matching one of these is automatically `DRIVER` |
| `loss_of_function_target` | `1` if any LoF SNV in this gene is `DRIVER` (e.g., TSGs) |
| `oncogenic_amplification` | `1` if a CNA amplification of this gene is `DRIVER` (e.g., CDK4, MDM2) |
| `lof_via_deletion` | `1` if a CNA homozygous deletion of this gene is `DRIVER` (e.g., CDKN2A, MTAP) |
| `notes` | Free-text rationale; quotes the case study or paper that motivates inclusion |

**The 21 genes** span the major RMS biology groups:

- Fusion oncoproteins: PAX3, PAX7, FOXO1
- Cell cycle / RB axis: CDK4, CDK6, CDKN2A, MTAP, MDM2
- RAS-MAPK: NRAS, KRAS, HRAS, NF1
- PI3K-AKT-mTOR: PTEN, PIK3CA
- RTKs: FGFR4
- Transcription / epigenetics: MYOD1, BRD4, BRD9, CHD4, PRMT5
- Stress response: TP53

Every change to this file is treated as an ADR-worthy scope decision. The targets list directly determines what events the pipeline can call as `DRIVER`, so adding or removing a gene is a science decision, not a code decision.

### 4.2 `depmap_rms_summary.tsv`: the dependency layer (PLACEHOLDER)

[assets/depmap_rms_summary.tsv](../assets/depmap_rms_summary.tsv). One row per gene with summary statistics across approximately 12 RMS cell lines. Schema:

| Column | What it does |
|---|---|
| `n_lines` | RMS lines surveyed |
| `mean_chronos_all` | Mean Chronos score across all RMS lines (lower is more essential; see [glossary](#a3-chronos-and-ceres-the-depmap-essentiality-scores)) |
| `mean_chronos_fp` | Subtype-stratified, fusion-positive RMS only |
| `mean_chronos_fn` | Subtype-stratified, fusion-negative RMS only |
| `pct_essential` | Fraction of lines with Chronos < -0.5 |
| `dependency_score` | Bounded 0-1 derived score (1 = strongest dependency, 0 = none) |
| `provenance` | Short note on why the values look the way they do |

This file is **PLACEHOLDER** in v0.4. The values are hand-curated to be consistent with published patterns:

- Dharia 2021 (Pediatric DepMap, Nature Genetics): RMS lineage essentialities
- Olanich 2015 (CDK4 in FP-RMS)
- Marjon 2016, Kryukov 2016, Mavrakis 2016 (PRMT5 synthetic lethality in MTAP-null)
- Gryder 2017-2020 (BRD4 dependency in FP-RMS via PAX3-FOXO1 super-enhancers)

Replacing this file with a real DepMap 24Q2+ pull is the highest-priority remaining work. The replacement requires the same schema. Phase 3 reads this file by column name, so a real pull with the right columns is a swap, not a code change.

### 4.3 `drug_target_map.tsv`: the curated drug catalog

[assets/drug_target_map.tsv](../assets/drug_target_map.tsv). Hand-curated from the v2 COG proposal drug list, the Molecular Targets Platform spot-checks, and FDA pediatric labeling. Schema:

| Column | What it does |
|---|---|
| `gene` | Target gene |
| `drug` | Drug name (lower-cased) |
| `mechanism` | Specific mechanism string (e.g., `CDK4/6_inhibitor`, `MEK_inhibitor`, `MTA_cooperative_PRMT5_inhibitor`) |
| `max_phase` | `approved` / `phase3` / `phase2` / `phase1` / `preclinical` |
| `pediatric_evidence` | `yes_approved` / `yes_trial` / `adult_only` / `none` |
| `notes` | Trial reference or rationale |

**The mechanism strings carry meaning.** The case-study scorecard asserts on substrings like `FGFR_inhibitor`, `MEK_inhibitor`, `CDK4/6_inhibitor`, `PRMT5_inhibitor`, `BET_inhibitor`. Curated entries are written to match those patterns. DGIdb-sourced rows use the more generic DGIdb labels (`inhibitor`, `antagonist`), which is why the curated map wins on (gene, drug) collisions in Phase 4.

### 4.4 `dgidb_drugs.tsv`: the live DGIdb cache (real, refreshable)

[assets/dgidb_drugs.tsv](../assets/dgidb_drugs.tsv). Auto-generated by [bin/fetch_dgidb.py](../bin/fetch_dgidb.py) from DGIdb's public GraphQL API. Same schema as `drug_target_map.tsv` so Phase 4 can union them.

The fetcher drops rows with no annotated mechanism type. Without this filter, drugs with mere literature co-mention (e.g., doxorubicin / MYOD1 from chemotherapy literature) flooded the report and pushed legitimate hits out of the top 10. The scorecard caught the regression on the first naive merge. After filtering, 1,381 raw interactions reduced to 198 mechanism-typed ones.

To refresh:

```bash
bin/fetch_dgidb.py
git add assets/dgidb_drugs.tsv
git commit -m "Refresh DGIdb cache"
```

---

## 5. The five phases, in depth

Each phase has the same shape: input file(s), what it does, output file, the algorithm or knowledge it depends on, and what's deferred to a later version.

### 5.1 Phase 1: variant + CNA + fusion annotation

**Code**: [bin/phase1_annotate.py](../bin/phase1_annotate.py), [modules/phase1_variants/main.nf](../modules/phase1_variants/main.nf).

**Inputs**:

- VCF (required): germline-filtered somatic SNV/indel calls. The toy fixtures put `GENE`, `CONSEQUENCE`, `EXPECTED`, and `NOTE` in the INFO column; a real pipeline will populate these from VEP + OncoKB.
- CNA TSV (optional): `sample_id, gene, event, copy_number, notes`. `event` is one of `amplification`, `homozygous_deletion`, `focal_gain`, `focal_loss`.
- Fusion TSV (optional): `sample_id, gene_5p, gene_3p, fusion_name, fusion_class, notes`.
- `assets/targets_kb.tsv`: the gene knowledge base.

**Output**: a unified per-event TSV with one row per molecular event. Column contract:

```
sample_id, event_id, event_type, chrom, pos, ref, alt,
gene, uniprot, role,
consequence, hgvsp_short, copy_number, fusion_partner,
is_target, call, reason, variant_score
```

`event_type` is `snv` / `cna` / `fusion`. `call` is `DRIVER` / `VUS` / `PASSENGER` / `OFF_TARGET`. `variant_score` is the first scoring component (1.0 for confident drivers, 0.4-0.6 for VUS, 0 for passengers).

**The classifier**: a small set of rules in `classify_snv`, `classify_cna`, `classify_fusion`:

- SNV LoF (stop_gained, frameshift, splice acceptor/donor) on a gene where `loss_of_function_target=1`: `DRIVER`, score 1.0.
- SNV missense matching a known hotspot from `targets_kb.tsv`: `DRIVER`, score 1.0.
- SNV protein-altering on a target gene but not a hotspot: `VUS`, score 0.4-0.6.
- SNV non-protein-altering (synonymous, intronic): `PASSENGER`, score 0.
- CNA amplification on a gene where `oncogenic_amplification=1`: `DRIVER`.
- CNA homozygous deletion on a gene where `lof_via_deletion=1`: `DRIVER`.
- CNA on a target gene at sub-threshold dosage: `VUS`.
- Fusion involving a gene where `role=fusion_partner`: `DRIVER` for both partners.
- Anything on a gene not in `targets_kb.tsv`: `OFF_TARGET`, dropped from downstream attention.

**The hotspot extractor** uses a regex tuned to the toy VCF's NOTE field convention (`R175H_canonical_hotspot_pathogenic` etc.). The regex is `(?<![A-Za-z0-9])([A-Z]\d+[A-Z*X])(?![A-Za-z0-9])` to handle underscores as delimiters. (This is a v0.x simplification. Real Phase 1 will rely on VEP's HGVSp output, not text mining.)

**Concept zoom**: see [A.1 Drivers vs passengers](#a1-driver-vs-passenger-variants), [A.2 Hotspots and loss of function](#a2-hotspot-residues-and-loss-of-function), [A.5 Copy number alterations](#a5-copy-number-alterations), [A.6 Fusion oncoproteins](#a6-fusion-oncoproteins).

**What's deferred to v0.5+**:

- Real VEP 111 + OncoKB / ClinVar / COSMIC annotation. The toy VCF supplies these as INFO fields for now.
- gnomAD-AF filtering for likely-germline variants.
- GISTIC2 / cnvkit / FACETS for real CNA calling from raw data.
- STAR-Fusion / Arriba for real fusion calling from RNA-seq.

### 5.2 Phase 2: structural reference attachment

**Code**: [bin/phase2_structure.py](../bin/phase2_structure.py), [modules/phase2_structure/main.nf](../modules/phase2_structure/main.nf).

**Input**: Phase 1 TSV.

**Output**: same rows plus four new columns:

```
alphafold_url, alphafold_pdb_url, structural_score, structural_reason
```

**What it does**: for each event on a target gene with a UniProt ID, attaches the AlphaFold DB entry URL (deterministic: `https://alphafold.ebi.ac.uk/entry/{uniprot}`) and the direct PDB download URL. These URLs are real and resolve to the actual reference structures.

**The structural rubric** (the second scoring component, weight 0.15):

- SNV at a known hotspot residue: 1.0. The reference structure resolves the residue and the literature already says the residue matters.
- SNV LoF in a TSG: 0.4. There is no folded mutant protein. The reference structure is still useful for context (e.g., where the truncation cuts off), but the structural-match concept does not apply.
- SNV protein-altering missense at a non-hotspot residue (VUS): 0.6. The reference is available, mutant prediction would tell you more but we have not done it yet.
- CNA amplification: 0.5. Gene-level event; the reference structure of the protein product is unaffected by dosage.
- CNA homozygous deletion: 0.4. No protein product; reference still useful for context.
- Fusion: 0.7. Both partner reference structures exist. Predicting the actual fusion junction structure (what the chimeric protein looks like) is the v0.5+ work.

**Concept zoom**: see [A.7 AlphaFold and structure prediction](#a7-alphafold-and-protein-structure-prediction), [A.8 Why structure helps drug ranking](#a8-why-protein-structure-helps-drug-ranking).

**What's deferred**:

- Boltz-1 (MIT) or Chai-1 (Chai Discovery) mutant-structure prediction for hotspot residues.
- AlphaFold-Multimer for fusion junction modelling (PAX3-FOXO1, PAX7-FOXO1).
- AutoDock Vina or Boltz ligand-aware mode for actual drug-protein docking.

These are GPU-bound and only worth doing on a shortlist of target-drug pairs, not the full corpus. v0.5 will add them on a configurable shortlist.

### 5.3 Phase 3: dependency integration

**Code**: [bin/phase3_dependency.py](../bin/phase3_dependency.py), [modules/phase3_dependency/main.nf](../modules/phase3_dependency/main.nf).

**Inputs**: Phase 2 TSV, `assets/depmap_rms_summary.tsv`, `--subtype` (FP, FN, UNKNOWN, or ALL).

**Output**: same rows plus seven new columns:

```
depmap_n_lines, depmap_chronos_all, depmap_chronos_fp, depmap_chronos_fn,
depmap_pct_essential, dependency_score, dependency_reason
```

**What it does**: joins each event's gene against the DepMap summary, picks the subtype-appropriate Chronos column, and computes a dependency_score in 0-1.

**The score computation**: starts from the precomputed `dependency_score` in the bundled summary, then adds a subtype-aware bonus when the requested subtype is more dependent than the all-RMS mean. Concretely:

```
chronos_used = depmap[subtype] (e.g., mean_chronos_fp if --subtype FP)
chronos_all  = depmap[ALL]
bonus = max(0, (chronos_all - chronos_used) * 0.3)
score = min(1.0, base_score + bonus)
```

This means a gene that is uniformly essential across RMS gets the base score; a gene that is preferentially essential in the requested subtype gets a bonus capped at 1.0.

Worked example. CDK4 in the bundled summary: `mean_chronos_all = -0.55`, `mean_chronos_fp = -0.80`, `base_score = 0.60`. For an FP-RMS run:
- `chronos_used = -0.80`, `chronos_all = -0.55`
- `bonus = max(0, (-0.55 - -0.80) * 0.3) = 0.075`
- `score = min(1.0, 0.60 + 0.075) = 0.675`

For an FN-RMS run, `chronos_used = -0.30`, bonus is 0 (FN is less dependent than ALL), `score = 0.60`. So CDK4 scores higher in FP than FN, which mirrors the published biology that CDK4 amplification is enriched in FP-RMS.

**Concept zoom**: see [A.3 Chronos and CERES](#a3-chronos-and-ceres-the-depmap-essentiality-scores), [A.4 CRISPR essentiality screens](#a4-genome-wide-crispr-essentiality-screens), [A.10 Subtype-stratified dependency](#a10-why-fp-vs-fn-rms-subtype-matters-for-dependency).

**What's deferred**:

- Real DepMap 24Q2+ pull. The bundled summary is hand-curated to match published patterns. The signed-URL flow for downloading from `storage.googleapis.com/depmap-external-downloads/...` requires authenticated portal API calls; the `depmap` Python package wraps this but adds a dependency.
- OpenPedCan expression z-scores. The pilot plan §4 Aim 2 Phase 3 calls for "expression z-score in RMS vs normal pediatric tissues" as a separate evidence dimension. This becomes the (currently deferred) expression-specificity component in Phase 5 scoring.
- Drug-sensitivity association from PPTC PDX models.

### 5.4 Phase 4: drug matching

**Code**: [bin/phase4_drugs.py](../bin/phase4_drugs.py), [modules/phase4_drugs/main.nf](../modules/phase4_drugs/main.nf).

**Inputs**: Phase 3 TSV, `assets/drug_target_map.tsv` (curated, primary), `assets/dgidb_drugs.tsv` (DGIdb cache, secondary).

**Output**: a long-format TSV with one row per (event, drug) pair, plus six new columns:

```
drug, drug_mechanism, drug_max_phase, drug_pediatric_evidence,
drug_evidence_score, drug_notes
```

**What it does**:

1. Loads both drug maps and unions them, keyed by `(gene, drug.lower())`. Curated wins on collisions because curated mechanism strings (`MEK_inhibitor`, `CDK4/6_inhibitor`, etc.) drive the case-study scorecard assertions.
2. For each Phase 3 row that is a `DRIVER` or `VUS`, looks up all drugs targeting that gene.
3. Applies conditional rules. v0.4 has one: `KRAS_G12C_inhibitor` (adagrasib, sotorasib) is suppressed unless the matched variant's `hgvsp_short` is exactly `G12C`.
4. Emits one row per applicable drug. Variants with no applicable drug (passengers, off-target, or genes with no mapped drug) get a single row with empty drug fields, so downstream variant counts stay honest.

**The drug_evidence_score** (the fifth scoring component, weight 0.20):

```
drug_evidence_score = phase_weight * pediatric_weight

phase_weight     = approved 1.0 | phase3 0.85 | phase2 0.70 | phase1 0.55 | preclinical 0.30
pediatric_weight = yes_approved 1.0 | yes_trial 0.85 | adult_only 0.65 | none 0.40
```

A single drug can contribute up to 1.0 to this component (FDA-approved with pediatric approval, like selumetinib for NF1 plexiform neurofibromas). A typical adult-only approved drug like erdafitinib scores 0.65. A phase-1 adult-only drug like MRTX1719 scores 0.36.

**Concept zoom**: see [A.11 DGIdb and drug-gene interaction databases](#a11-dgidb-and-drug-gene-interaction-databases), [A.12 Pediatric drug evidence levels](#a12-pediatric-drug-evidence-levels-why-yes_approved--adult_only).

**What's deferred**:

- OpenTargets and the Molecular Targets Platform (MTP) for additional target-drug evidence and for pediatric-specific weighting.
- ClinicalTrials.gov for live trial-stage status.
- CMap / LINCS L1000 for signature-based drug matching (project the tumor's transcriptional profile onto perturbation signatures, retrieve compounds that reverse it).
- Conditional rules beyond G12C: MDM2 inhibitors should require TP53 wild-type, BET inhibitors should weight up on PAX3-FOXO1 fusion presence, etc. These will move into a `condition` column on the drug map rather than being hard-coded.

### 5.5 Phase 5: confidence scoring and report

**Code**: [bin/phase5_score.py](../bin/phase5_score.py), [modules/phase5_scoring/main.nf](../modules/phase5_scoring/main.nf).

**Inputs**: Phase 4 long-format TSV, the original VCF (used for the report's input fingerprint), `--sample-id`, `--subtype`, `--pipeline-version`.

**Outputs**:

- `<sample>.phase5.tsv`: the long-format scored TSV with `confidence` and per-component columns added.
- `<sample>.report.md`: human-readable Markdown report ranked by confidence.

**The confidence formula** (the heart of the pipeline):

```
confidence = 0.25 * variant_score
           + 0.15 * structural_score
           + 0.25 * dependency_score
           + 0.15 * expression_score    [DEFERRED in v0.4, contributes 0]
           + 0.20 * drug_evidence_score
```

Weights come from the pilot plan §4 Aim 2 Phase 5. They are heuristic, fixed by prior, and documented as such in the report itself. The pilot plan §11 Risk 5 explicitly accepts that the first version's formula is heuristic; the scientific value is that every component is individually defensible and visible.

**Why expression is deferred**: Phase 3 v0.4 only delivers a dependency score, not an expression-specificity score. The pilot plan separates these because they answer different questions: dependency = "does the cell line need this gene to live"; expression specificity = "is this gene over-expressed in RMS vs normal pediatric tissues". v0.5 will add OpenPedCan expression integration to fill the 0.15 weight.

**Why a heuristic over a learned model**: per the pilot plan §11 Risk 5, "the first version is heuristic with fixed weights ... a second version, calibrated against MCI data and held-out validation tumors, can be learned rather than hand-tuned, but only after we know the pipeline works." Learning weights without holdout data is overfitting masquerading as calibration.

**The report structure**:

1. Header: sample, subtype, pipeline version, run timestamp, input VCF SHA, event count breakdown.
2. v0.x disclaimer pointing out the placeholders.
3. Top-10 ranked target-drug pairs.
4. Per-event detail (sorted by best confidence): event metadata, all four component scores, dependency context line, AlphaFold link, then a per-drug ranking table for that event.
5. Methodology section with the weights printed in plain text.
6. Footer: "engineering validation output, not medical advice".

The per-event detail block is what makes the pipeline auditable. A clinician can land on any row of the top-10 and see exactly why it ranked there, including the failure modes (which component dragged it down or pulled it up).

**Concept zoom**: see [A.13 Why a transparent scoring formula matters](#a13-why-a-transparent-scoring-formula-matters), [A.14 Score calibration vs ranking](#a14-score-calibration-vs-ranking).

---

## 6. End-to-end: what happens when you press go

### 6.1 The default invocation

```bash
nextflow run main.nf -profile laptop
```

This runs against the default fixture (TOY_TUMOR, FN-RMS, SNV-only) using `nextflow.config` defaults. The `-profile laptop` selects local execution with no Docker (Docker support is wired but not enabled until containers exist).

### 6.2 What happens, in order

1. Nextflow parses `nextflow.config`, loads param defaults, and includes `conf/laptop.config`.
2. `main.nf` builds five channels: the input VCF, the CNA TSV (default `assets/empty.cna.tsv`), the fusion TSV (default `assets/empty.fusion.tsv`), the targets KB, the DepMap summary, the curated drug map, and the DGIdb cache.
3. `PHASE1_ANNOTATE` is invoked with the sample inputs and the targets KB. Inside the process, `phase1_annotate.py` runs and writes `TOY_TUMOR.phase1.tsv` to the work directory. Nextflow stages this into `results/phase1/`.
4. `PHASE2_STRUCTURE` consumes the Phase 1 output and adds AlphaFold + structural columns. Output published to `results/phase2/`.
5. `PHASE3_DEPENDENCY` consumes the Phase 2 output plus the DepMap summary and the subtype param. Output to `results/phase3/`.
6. `PHASE4_DRUGS` consumes the Phase 3 output plus both drug maps. Output to `results/phase4/`.
7. `PHASE5_SCORE` consumes the Phase 4 output plus the original VCF (for the input fingerprint) and writes both the scored TSV and the markdown report to `results/phase5/`.
8. Nextflow prints `report ready: TOY_TUMOR -> ...` and exits.
9. `results/pipeline_info/` gets the trace, timeline, report, and DAG artifacts.

Total wall time on a laptop: roughly 5-10 seconds. No GPU, no large downloads.

### 6.3 Running a different fixture

For the FP-RMS toy patient with CDK4 amp + PAX3-FOXO1 fusion:

```bash
nextflow run main.nf -profile laptop \
    --input  tests/data/toy_fp_cdk4amp.vcf \
    --cna    tests/data/toy_fp_cdk4amp.cna.tsv \
    --fusion tests/data/toy_fp_cdk4amp.fusion.tsv \
    --sample_id TOY_FP_CDK4 \
    --subtype FP
```

For all three toy patients with a cross-sample summary plus the scorecard:

```bash
bin/run_all_toys.sh
open results/multisample_summary.md
open results/scorecard.md
```

### 6.4 The output structure per sample

```
results/
├── phase1/<sample>.phase1.tsv          per-event annotation
├── phase2/<sample>.phase2.tsv          + structure
├── phase3/<sample>.phase3.tsv          + dependency
├── phase4/<sample>.phase4.tsv          + drugs (long format)
├── phase5/
│   ├── <sample>.phase5.tsv             + confidence + per-component scores
│   └── <sample>.report.md              the human-readable report
└── pipeline_info/
    ├── trace.txt
    ├── report.html
    ├── timeline.html
    └── dag.svg
```

Every TSV has a stable column contract. Every TSV is a strict superset of its input's columns. This is what makes substituting any phase a column-conformant change rather than a refactor.

---

## 7. How we know it works: the case-study scorecard

### 7.1 Why a scorecard

The pilot plan §3 defines four scientific case studies. Each one is a known piece of RMS biology with a known correct answer. The pipeline is trustworthy on new data only if it produces the correct answer on these known cases.

The scorecard turns that judgment into a machine-checkable test. It runs the pipeline on each fixture and asserts that the expected drug class shows up at the expected rank. CI runs the scorecard on every push and blocks merge if any assertion regresses.

### 7.2 The four case studies plus the passenger sanity check

Encoded as `[[case]]` entries in [tests/cases.toml](../tests/cases.toml):

| ID | Pilot # | Tumor profile | Asserts |
|---|---|---|---|
| `case3_fgfr4` | 3 | FGFR4 V550L hotspot in TOY_TUMOR | `FGFR_inhibitor` in top 5 |
| `case4_ras_mek` | 4 | NRAS Q61K + KRAS G12V in TOY_TUMOR | `MEK_inhibitor` in top 2 |
| `toy_tumor_passengers` | general | 3 passengers in TOY_TUMOR | All passengers below 0.10 confidence |
| `case2_cdk4_amp` | 2 | CDK4 amp + PAX3-FOXO1 in TOY_FP_CDK4 | `CDK4/6_inhibitor` in top 3, `CDK4` gene in top 3, `BET_inhibitor` in top 6 |
| `case1_mtap_prmt5` | 1 | MTAP+CDKN2A 9p21 co-deletion in TOY_MTAP_NULL | `PRMT5_inhibitor` in top 7, `CDK4/6_inhibitor` in top 3 |

8 assertions, 5 cases. All currently PASS at v0.4.

### 7.3 The three assertion kinds

Implemented in [bin/check_case_studies.py](../bin/check_case_studies.py):

- **`mechanism_in_top_n`**: drug_mechanism column in top-N rows must contain the given substring. Substring match means `FGFR_inhibitor` matches both `pan_FGFR_inhibitor` and `FGFR4_selective` (which is intentional: both are "FGFR inhibitors" for case-study purposes).
- **`gene_event_in_top_n`**: a specific gene must appear in the top-N rows. Used for case 2 to assert that CDK4 itself (not just CDK4/6 inhibitors against CDKN2A loss) shows up.
- **`passenger_below`**: every PASSENGER-call row must have confidence below the threshold. Catches regressions where the passenger filter breaks.

### 7.4 What "PASS" means and what it doesn't

PASS means the pipeline produces the right top-N drug for known biology. It does NOT mean:

- The dependency scores are correct in absolute terms (Phase 3 is still a placeholder).
- The drug list is exhaustive (curated map plus DGIdb is not the whole pharmacopeia).
- The pipeline will discover new biology (the pilot plan §7 is explicit: getting known biology right just establishes a floor).

PASS is necessary but not sufficient for trusting the pipeline on real MCI data.

### 7.5 The CI workflow

[.github/workflows/ci.yml](../.github/workflows/ci.yml) has three jobs:

1. **`case-study-scorecard`**: Python-only, fastest, the primary correctness gate. Uses Python 3.11 (for `tomllib`) and runs `python3 bin/check_case_studies.py --quiet`. Uploads `results/scorecard.{md,json}` as artifacts.
2. **`nextflow-end-to-end`**: installs Java 17 + Nextflow, runs the pipeline on TOY_TUMOR (default fixture) and TOY_FP_CDK4 (CNA + fusion fixture), verifies the expected report files exist.
3. **`markdown-lint`**: greps every committed `*.md` for em dashes and fails if any are found. Project style rule.

---

## 8. What's real, what's placeholder, what's deferred

| Phase | v0.4 status | What's real | What's PLACEHOLDER | What's deferred entirely |
|---|---|---|---|---|
| 1 | Curated KB join | The hotspot list, LoF flags, amp/del flags in `targets_kb.tsv` are real biology | Toy VCF supplies `GENE`/`CONSEQUENCE` in INFO; real VEP not invoked | gnomAD-AF germline filter, CNV calling from BAMs, fusion calling from RNA-seq |
| 2 | Reference attachment | AlphaFold DB URLs are real and resolve to actual structures | The structural_score is heuristic, not docking-derived | Boltz-1 / Chai-1 mutant prediction, AlphaFold-Multimer for fusion junctions, AutoDock Vina docking |
| 3 | Chronos summary join | Schema and join logic | Per-gene Chronos values are hand-curated, not from real DepMap | Real DepMap pull, OpenPedCan expression, PPTC drug response |
| 4 | Curated map UNION DGIdb | DGIdb GraphQL pull is live; curated mechanism strings are real | Curated map is not exhaustive | OpenTargets, MTP, ClinicalTrials.gov, CMap/LINCS L1000 |
| 5 | Weighted-sum scoring | The formula and the weights are the published ones from the pilot plan | The expression-specificity component is held at 0 | Learned weights from MCI data |

### 8.1 Why this much placeholder is OK for v0.x

The pilot plan §4 Aim 1 calls for a runnable end-to-end pipeline before chasing data fidelity. Building all five phases first, with bundled curated data, lets us:

- Validate the architecture (the column contracts, the Nextflow wiring, the report shape) cheaply.
- Catch integration bugs when each phase is small.
- Run the case-study scorecard from day one as a regression-prevention tool.
- Onboard new collaborators in under an hour.

When we replace placeholders with real data, the scorecard is what tells us whether the swap broke anything scientific. It already caught one regression during DGIdb integration (interaction-only rows flooded the report with literature noise; the scorecard caught it on the first run).

---

## 9. Where the pipeline goes from here

In rough priority order:

### 9.1 Real DepMap data (highest priority)

Replaces the biggest remaining placeholder. DepMap publishes per-release CRISPR Chronos scores via a signed-URL flow served by the portal API. Two paths:

- Use the `depmap` Python package (adds a dependency).
- Reverse-engineer the signed-URL flow from the portal and use plain Python.

The output is the same `assets/depmap_rms_summary.tsv` schema, just with real values. A `bin/fetch_depmap.py` script analogous to `bin/fetch_dgidb.py` would refresh it on demand.

### 9.2 Live OpenTargets and ClinicalTrials.gov for Phase 4

Same pattern as DGIdb: a fetcher writes a cached TSV, Phase 4 unions all sources with curated winning. OpenTargets adds disease-association evidence and pediatric-trial status. ClinicalTrials.gov adds live trial-stage information.

### 9.3 Containerize all five phases

The Dockerfiles in `containers/` are Ubuntu placeholders. Filling them in with `python:3.11-slim` (sufficient for v0.4, which has no real bioinformatics tools yet) gets `docker.enabled = true` working in the laptop profile. Real Phase 1 will need a VEP container, which is a much bigger build.

### 9.4 Real VEP + OncoKB for Phase 1

The toy VCFs supply `GENE`/`CONSEQUENCE` in INFO so the v0.x classifier can ignore the annotation problem. Real Phase 1 will run VEP 111 against GRCh38 + GENCODE v44 with ClinVar, COSMIC Cancer Gene Census, OncoKB, AlphaMissense, REVEL, and CADD. This is heavy tooling (VEP cache is ~50GB) and only fits a containerized phase, not the laptop profile.

### 9.5 Aim 1 corpus assembly (the first real-data run)

Pull processed MAF files for a few TARGET-RT samples from the GDC public tier, run them through the pipeline. First time the pipeline sees a real RMS tumor. Validates that the column contract holds when inputs come from outside our control.

### 9.6 v0.5 structural modelling on a shortlist

Phase 2 v0.5 will integrate Boltz-1 or Chai-1 to predict mutant structures for the top-N candidates per tumor (configurable), and AutoDock Vina or Boltz ligand-aware mode for actual docking. GPU-bound, so it gates on the `gpu` profile being live.

### 9.7 The MCI handoff

The pilot plan §12 documents what changes when MCI access arrives:

- Input file paths via `config/mci.yaml` (no code change).
- Reference genome: MCI uses GRCh38, pilot does too. No build translation.
- Sample metadata: MCI uses CCDI data model. Phase 1 ingestion needs a CCDI-to-internal-schema adapter (exercised against public CCDI submissions during pilot).
- QC thresholds re-tuned on the first 10-20 MCI samples.
- HIPAA-compliant compute (BCM RICC or institutional cloud).
- IRB protocol covering secondary analysis of de-identified dbGaP data.

The pilot's validation artifacts (benchmark tables, case study outputs, confidence calibration) become the regression tests that the MCI run must not break.

---

## Appendix A: glossary of biological and computational concepts

Each entry stands alone. Voice-chat Claude can jump to any of these when a term comes up.

### A.1 Driver vs passenger variants

A **driver** mutation gives the tumor a fitness advantage. A **passenger** mutation rides along but does not contribute to oncogenesis. The pipeline's first job is to separate the two.

The classification rules are mostly biological convention:

- A missense at a recurrent hotspot residue (e.g., NRAS Q61K, KRAS G12V, BRAF V600E) is almost always a driver because the residue is structurally critical.
- A truncating mutation (stop_gained, frameshift) in a tumor suppressor gene is a driver because it ablates function.
- A synonymous variant or an intronic variant outside splice regions is almost always a passenger.

The hard cases are non-recurrent missense mutations in target genes, which we call VUS (variant of unknown significance). The pipeline scores them at 0.4-0.6 to reflect uncertainty.

### A.2 Hotspot residues and loss of function

A **hotspot** is a specific amino acid residue where a mutation recurrently activates an oncogene. Famous examples:

- TP53 R175H, R248Q, R273H: structural hotspots; the mutant protein misfolds and loses TSG function (and gains some neomorphic activity).
- NRAS / KRAS / HRAS Q61, G12, G13: GTPase-cycle hotspots; mutation locks the protein in the active GTP-bound state.
- BRAF V600E: kinase-domain hotspot; constitutively active monomer.
- MYOD1 L122R: a neomorphic hotspot specific to a high-risk FN-RMS subtype (Shern 2014).
- FGFR4 V550L, N535K: kinase-domain hotspots in roughly 10% of RMS (Shern 2014).

**Loss of function (LoF)** is the opposite mode: any mutation that destroys the protein's function counts. Stop_gained, frameshift, and splice-disrupting mutations are LoF by default. LoF mutations matter most in tumor suppressor genes (TP53, NF1, CDKN2A, PTEN), where losing function unleashes growth.

The two modes are mutually exclusive at the score level: a hotspot residue and a stop_gained on the same gene get the same DRIVER call, but the reasoning is different (gain-of-function vs loss-of-function).

### A.3 Chronos and CERES, the DepMap essentiality scores

The DepMap project does genome-wide CRISPR knockouts in over 1,000 cancer cell lines, then measures how each cell line responds to losing each gene. The score per (gene, cell line) reflects how essential the gene is for that line's survival.

- **CERES** was the original score, normalized to 0 = essentiality of a known nonessential gene set, -1 = essentiality of a known pan-essential gene set.
- **Chronos** is the current generation, similar interpretation but with better statistical modeling of the screen dynamics. Lower is more essential.

Convention:
- Chronos > 0: the gene is dispensable for that line.
- Chronos near -0.5: the gene is essential.
- Chronos < -1: the gene is strongly essential.

The Pediatric DepMap (Dharia 2021, Nature Genetics) extended these screens to about a dozen RMS lines covering both fusion-positive and fusion-negative lineages. Examples: BRD4 strongly essential pan-RMS, CDK4 essential preferentially in FP-RMS (where 12q amplification is recurrent), PAX3 essential in FP-RMS via the PAX3-FOXO1 fusion.

### A.4 Genome-wide CRISPR essentiality screens

The screen design: a library of single guide RNAs (sgRNAs), four to six per gene, covering all roughly 18,000 protein-coding genes, is delivered to a pool of cells via lentivirus at low multiplicity (one sgRNA per cell). Cas9 cuts the targeted locus, the cell's repair machinery introduces frameshift indels, the gene is functionally knocked out in that cell. The pool grows for 14-21 days. sgRNA abundance is sequenced before and after. Guides targeting essential genes drop out (the cells died); guides targeting nonessential genes stay constant. The dropout per gene, normalized for sgRNA efficiency and copy-number bias, is the essentiality score.

Why this matters for drug discovery: a gene that scores essential in a cell line is a candidate target for inhibition in a tumor of that lineage. The dependency score in our Phase 3 is a proxy for "would knocking out this gene therapeutically work".

### A.5 Copy number alterations

DNA copy number is normally 2 (one chromosome per parent). Cancer cells routinely have segments amplified (extra copies, sometimes 5-50x) or deleted (one or both copies lost). Two clinically actionable patterns in RMS:

- **CDK4 amplification at 12q13-14**: about 10% of FP-RMS. The extra CDK4 copies overdrive the cell cycle, creating a CDK4/6 inhibitor opportunity (palbociclib, ribociclib, abemaciclib). Often co-amplifies MDM2.
- **CDKN2A homozygous deletion at 9p21**: about 15-20% of FN-RMS. Loses the p16 brake on CDK4/6, so paradoxically also creates a CDK4/6 inhibitor opportunity. Co-deletes MTAP (the next gene over on 9p21).

The MTAP co-deletion creates a synthetic lethality opportunity. See [A.9](#a9-synthetic-lethality-and-the-mtapprmt5-example).

### A.6 Fusion oncoproteins

A chromosomal translocation can fuse two genes into a chimeric protein with new function. The classic RMS example is t(2;13)(q35;q14), which fuses PAX3 (a transcription factor) to FOXO1 (another transcription factor). The chimeric PAX3-FOXO1 protein binds PAX3's DNA targets but with FOXO1's activation domain, driving an aberrant transcriptional program. PAX3-FOXO1 defines the FP-RMS subtype and is the founding genetic event in those tumors. About 10-15% of FP-RMS instead carries PAX7-FOXO1 from t(1;13)(p36;q14).

Targeting fusion oncoproteins directly is hard (no enzymatic activity to inhibit). Indirect targeting works via:

- **BET inhibitors** (birabresib, molibresib, mivebresib): BRD4 co-occupies 95% of PAX3-FOXO1's super-enhancers (Gryder 2017). BET inhibition disrupts that super-enhancer-driven transcription.
- **CDK7 inhibitors**: similar logic, hits the basal transcriptional machinery the fusion depends on.
- **Targeted protein degraders (BRD9 degraders CFT-8634, FHD-609)**: ablate the BAF complex member BRD9 that the fusion-bound chromatin needs.

### A.7 AlphaFold and protein structure prediction

AlphaFold 2 (DeepMind, 2021) predicts 3D protein structures from amino acid sequence with accuracy comparable to experimental crystallography for many proteins. The AlphaFold Database (alphafold.ebi.ac.uk) hosts predictions for over 200 million proteins, including the entire human proteome. Each entry has a confidence score per residue (pLDDT) and overall.

Why we use it in Phase 2: for any target gene, we can pull the reference structure for free without crystallography. The URL is deterministic (`https://alphafold.ebi.ac.uk/entry/{uniprot}`), so Phase 2 just attaches the URL plus a file pointer to the .pdb model.

What AlphaFold does NOT do well (yet):

- Multi-domain proteins with flexible linkers (the predicted poses are often arbitrary).
- Conformational changes (it predicts one state, not the ensemble).
- Mutant structures (the database stores the wild-type sequence; predicting a mutant requires re-running the model).
- Fusion junctions (chimeric proteins are not in the database).

Boltz-1 (MIT) and Chai-1 (Chai Discovery) are next-gen open-source predictors that handle on-demand mutant prediction and ligand-aware folding. Phase 2 v0.5 will use them on a shortlist.

### A.8 Why protein structure helps drug ranking

Drugs bind specific pockets on proteins. If a hotspot mutation sits at the binding pocket of a known drug, the mutation likely changes the drug's affinity (could go either way: better binding or escape). If the mutation is far from the pocket, the drug should still work.

In v0.4 we use structure as a confidence signal, not a docking score: the existence of a known reference structure for a hotspot residue means the residue is well-characterized, the literature exists, and the structural community has thought about it. That justifies the 1.0 structural_score for hotspots.

In v0.5 with real docking, structure becomes quantitative: predicted binding affinity for each (mutant, drug) pair, used to filter or rank.

### A.9 Synthetic lethality and the MTAP/PRMT5 example

**Synthetic lethality**: gene A and gene B are synthetically lethal if losing either alone is tolerated but losing both kills the cell. The therapeutic angle: if a tumor has lost gene A, a drug that inhibits gene B will kill the tumor selectively (normal cells still have A and survive).

The textbook example in our case: MTAP and PRMT5.

- MTAP is the methylthioadenosine phosphorylase. Normally it cleaves MTA (a byproduct of polyamine synthesis) into adenine and methylthioribose-1-phosphate.
- When MTAP is lost (typically via 9p21 co-deletion with CDKN2A), MTA accumulates in the cell.
- Accumulated MTA partially inhibits PRMT5 (a methyltransferase) by competing with its cofactor SAM.
- Cells with MTAP loss therefore have reduced baseline PRMT5 activity but still need PRMT5 to live.
- A drug that further inhibits PRMT5 in this MTAP-low context tips the balance to lethal. Normal cells (MTAP+, full PRMT5 activity) survive the same drug dose.
- MRTX1719 (Mirati) is a PRMT5 inhibitor specifically designed to be MTA-cooperative: it preferentially inhibits the PRMT5-MTA complex, achieving a wider therapeutic window in MTAP-null tumors.

This is exactly what pilot case study 1 tests: TOY_MTAP_NULL has CDKN2A + MTAP homozygous co-deletion, and the pipeline puts MRTX1719 in the top 7 alongside the CDK4/6 inhibitors that come from the CDKN2A loss.

### A.10 Why FP vs FN RMS subtype matters for dependency

FP-RMS and FN-RMS are biologically different diseases despite sharing histology:

- FP-RMS is fusion-driven (PAX3-FOXO1 or PAX7-FOXO1). The fusion creates a transcriptional addiction: the cell depends on continued PAX3-FOXO1 expression and the super-enhancer landscape it builds. CDK4 amplification is enriched. BRD4 is highly essential. PAX3 itself is essential (because the fusion locks it in).
- FN-RMS is mutation-driven (RAS-MAPK pathway, MYOD1 L122R, FGFR4 hotspots). No defining fusion. CDKN2A loss is more common. RAS hotspot mutations are more common.

The implication for dependency scoring: a gene that is essential in FP-RMS may be uninteresting in FN-RMS, and vice versa. Our `mean_chronos_fp` and `mean_chronos_fn` columns let Phase 3 weight the dependency by the requested subtype. The `--subtype` flag plumbs through to Phase 3 from `nextflow.config`.

When subtype is unknown (e.g., a real tumor we haven't classified yet), Phase 3 falls back to `mean_chronos_all`. v0.5 will add an automated subtype detector that reads the fusion input and sets subtype to FP if a driver fusion is present.

### A.11 DGIdb and drug-gene interaction databases

The Drug-Gene Interaction Database (DGIdb, dgidb.org) is a meta-aggregator. It pulls drug-gene relationships from approximately 30 source databases (CIViC, OncoKB, ChEMBL, GuideToPharmacology, DrugBank, NCI, MyCancerGenome, ClinicalTrials registry mentions, and many more) and harmonizes them into a single GraphQL API.

For each (drug, gene) pair, DGIdb gives:

- The drug name, RxNorm or NCIt or DrugBank concept ID, and approval status.
- One or more interaction types (inhibitor, antagonist, agonist, modulator, blocker, etc.) when annotated. Many interactions are not mechanistically annotated and are tagged only as "interaction".
- The list of source databases that contributed evidence.
- PubMed IDs for primary literature.

What DGIdb is good at: breadth. For a target like CDK4, DGIdb returns 50+ compounds, many you would not have curated by hand.

What DGIdb is bad at: pediatric-specific evidence (it does not flag pediatric trials). Mechanistic specificity (the labels are coarse, and uncategorized interactions are noise unless filtered).

Our Phase 4 strategy: union DGIdb with a hand-curated map. Curated wins on collisions because curated mechanism strings drive the case-study assertions. DGIdb adds breadth in the per-event detail tables. Mechanism-unannotated DGIdb rows are dropped to keep the report clean.

### A.12 Pediatric drug evidence levels: why yes_approved > adult_only

The drug map has a `pediatric_evidence` column with four values:

- **yes_approved**: FDA-approved in a pediatric indication. Selumetinib for NF1 plexiform neurofibromas is the canonical example. These score 1.0 on the pediatric-evidence weight.
- **yes_trial**: actively in pediatric clinical trials, not yet approved. COG palbociclib trial in sarcoma fits here. Score 0.85.
- **adult_only**: approved or in trials in adults but not yet pediatric. Erdafitinib (FGFR inhibitor approved for adult urothelial cancer) is here. Score 0.65.
- **none**: no clinical use in humans yet. Tool compounds, preclinical leads. Score 0.40.

Why this matters for ranking: a drug that has crossed the pediatric trial threshold is much closer to clinically useful for an RMS patient than an adult-only drug. The 0.65 vs 1.0 weight separation is intentional. Translating an adult-only nominee to a pediatric trial concept requires bridging adult PK/safety data, dose-finding work, and IRB approval, which can take years.

### A.13 Why a transparent scoring formula matters

Two reasons.

First, **clinician trust**. A weighted-sum confidence score with five named components, each visible in the per-event detail, is something a pediatric oncologist can read and either accept or push back on. A black-box neural network output of "0.73" is not. The pipeline's eventual customer is the COG STS committee; they need to be able to defend the recommendation in a tumor board.

Second, **pilot plan §11 Risk 6 discipline**: "If a case study returns wrong top drug class, find the phase-level bug, do not widen the confidence weights to rescue the expected answer." A transparent formula makes that bug-finding tractable. When a case fails, you can read the per-event detail and see which component is the culprit (the dependency score is too low because the cell-line panel lacks the relevant subtype, or the drug evidence weight is too low because pediatric data is missing, etc.). With a learned weight matrix you cannot diagnose this; you can only retrain.

The pilot plan is explicit that v1 may be a learned model, but only AFTER we have held-out validation tumors to calibrate against. That requires MCI data, which we do not have yet. For now, fixed weights from the pilot plan are correct.

### A.14 Score calibration vs ranking

The pilot plan acceptance criterion is "expected drug class in top-3 (or top-N)". This is a **ranking** criterion, not a **calibration** criterion. The pipeline does not need to produce well-calibrated absolute confidence scores; it needs to produce the correct rank order.

This matters because the deferred expression component (currently contributing 0) means absolute confidence scores cap at 0.85 instead of 1.00. That is a calibration shift, not a ranking shift, because it affects every drug equally. The case-study scorecard would still pass even if we renormalized the weights, because the relative ranking is preserved.

When we add real expression data in v0.5, absolute scores will go up, but the ranking should not flip. If it does, that signals a real biological tension (the expression evidence disagrees with the dependency evidence) and warrants investigation.

---

## Appendix B: file-by-file reference

Quick lookup. Voice-chat Claude can use this to point you at a specific file when discussing a topic.

### B.1 Top-level

- [README.md](../README.md): project overview, status, quickstart, eyeball-tests for case studies.
- [main.nf](../main.nf): top-level Nextflow workflow. Builds channels, calls each PHASE process in sequence, prints the final report path.
- [nextflow.config](../nextflow.config): all param defaults including input fixtures, subtype, asset paths. CLI `--<param>` overrides everything here.
- [LICENSES.md](../LICENSES.md): per-data-source license table.

### B.2 Configuration

- [conf/base.config](../conf/base.config): shared resource defaults, error retry strategy, label-based resource scaling (cpu_small, cpu_medium, cpu_large, gpu).
- [conf/laptop.config](../conf/laptop.config): local CPU profile. Docker disabled until containers exist.
- [conf/slurm.config](../conf/slurm.config): institutional HPC profile.
- [conf/aws.config](../conf/aws.config): AWS Batch profile.
- [conf/gpu.config](../conf/gpu.config): GPU profile for future Boltz/Chai work.
- [conf/targets.yaml](../conf/targets.yaml): the original 21-gene target list with per-gene rationale. The header warns that every change is ADR-worthy.

### B.3 Nextflow modules (thin wrappers)

- [modules/phase1_variants/main.nf](../modules/phase1_variants/main.nf): wraps `phase1_annotate.py`, takes (sample_id, vcf, cna, fusion) tuple plus targets_kb path.
- [modules/phase2_structure/main.nf](../modules/phase2_structure/main.nf): wraps `phase2_structure.py`.
- [modules/phase3_dependency/main.nf](../modules/phase3_dependency/main.nf): wraps `phase3_dependency.py`, threads subtype param.
- [modules/phase4_drugs/main.nf](../modules/phase4_drugs/main.nf): wraps `phase4_drugs.py`, accepts both curated and extra drug maps.
- [modules/phase5_scoring/main.nf](../modules/phase5_scoring/main.nf): wraps `phase5_score.py`, threads pipeline_version.

### B.4 Python implementations

- [bin/phase1_annotate.py](../bin/phase1_annotate.py): the SNV / CNA / fusion classifier. Heart of the rule logic.
- [bin/phase2_structure.py](../bin/phase2_structure.py): AlphaFold attachment + heuristic structural_score.
- [bin/phase3_dependency.py](../bin/phase3_dependency.py): DepMap join + subtype-aware score.
- [bin/phase4_drugs.py](../bin/phase4_drugs.py): drug map union + conditional rules + drug_evidence_score.
- [bin/phase5_score.py](../bin/phase5_score.py): weighted-sum confidence + Markdown report renderer.
- [bin/fetch_dgidb.py](../bin/fetch_dgidb.py): DGIdb GraphQL fetcher. Run manually to refresh the cache.
- [bin/check_case_studies.py](../bin/check_case_studies.py): case-study scorecard runner. Used by CI.
- [bin/run_all_toys.sh](../bin/run_all_toys.sh): runs all three toy patients via Nextflow plus the scorecard.

### B.5 Reference data

- [assets/targets_kb.tsv](../assets/targets_kb.tsv): the 21-gene knowledge base. Hotspots, LoF flags, amp/del flags, UniProt IDs.
- [assets/depmap_rms_summary.tsv](../assets/depmap_rms_summary.tsv): PLACEHOLDER Chronos summary. Highest-priority real-data swap.
- [assets/drug_target_map.tsv](../assets/drug_target_map.tsv): curated drug-target map with mechanism strings, max_phase, pediatric_evidence.
- [assets/dgidb_drugs.tsv](../assets/dgidb_drugs.tsv): auto-generated DGIdb cache (mechanism-typed only).
- [assets/empty.cna.tsv](../assets/empty.cna.tsv), [assets/empty.fusion.tsv](../assets/empty.fusion.tsv): default empty inputs for samples without CNA / fusion data.
- [assets/README.md](../assets/README.md): provenance and replacement schema documentation.

### B.6 Test fixtures and acceptance criteria

- [tests/cases.toml](../tests/cases.toml): the source of truth for case-study acceptance assertions.
- [tests/data/toy_patient.vcf](../tests/data/toy_patient.vcf): TOY_TUMOR fixture, FN-RMS, 7 driver SNVs + 3 passengers.
- [tests/data/toy_fp_cdk4amp.{vcf,cna.tsv,fusion.tsv}](../tests/data/): TOY_FP_CDK4 fixture, FP-RMS with PAX3-FOXO1 + CDK4/MDM2 amp + TP53 R175H.
- [tests/data/toy_mtap_null.{vcf,cna.tsv,fusion.tsv}](../tests/data/): TOY_MTAP_NULL fixture, FN-RMS with CDKN2A + MTAP 9p21 co-deletion.
- [tests/data/README.md](../tests/data/README.md): describes each fixture and how to invoke the pipeline against it.

### B.7 Documentation

- [docs/walkthrough.md](walkthrough.md): this file.
- [docs/architecture.md](architecture.md): short note pointing at the pilot plan.
- [docs/data_sources.md](data_sources.md): placeholder; populated at Aim 1 close with retrieval scripts and corpus manifest schemas.
- [docs/MCI_TRANSITION.md](MCI_TRANSITION.md): placeholder; the configuration delta for the MCI handoff.
- [docs/adr/0001-record-architecture-decisions.md](adr/0001-record-architecture-decisions.md): the ADR convention itself.
- [docs/adr/0002-nextflow-orchestrator.md](adr/0002-nextflow-orchestrator.md): why Nextflow over Snakemake.

### B.8 CI

- [.github/workflows/ci.yml](../.github/workflows/ci.yml): three jobs (case-study scorecard, Nextflow end-to-end, em-dash markdown lint).

### B.9 The plans (gitignored, in main repo dir)

- `plans/pilot_ccdi_mci/rms_translational_pilot_project.md`: the canonical build plan. 4 Aims, 5 phases, 4 case studies, 12-month timeline, explicit success criteria.
- `plans/pilot_ccdi_mci/ccdi_mci_data_overview.md`: the CCDI / MCI / TARGET / DepMap / MTP data landscape with access-tier matrix.
- `plans/v2_rms_in_silico_pipeline_proposal.docx`: the parent COG concept proposal.
- `plans/v2_expanded_rms_in_silico_pipeline.docx`: the expanded version of v2.
- `plans/AI_Translational_Therapeutics_Platform_v2.docx`: the umbrella program.
