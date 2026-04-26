# ADR 0003: Pluggable Annotators with API-First Defaults

Date: 2026-04-26

## Status

Accepted.

## Context

The pipeline will be applied to many cohorts and eventually to other diseases. Variant annotation is the first phase that wants more than one backend: toy fixtures use hand-curated `GENE` / `CONSEQUENCE` baked into VCF INFO; real cohorts need real annotation; future work will overlay OncoKB driver classification, ClinVar pathogenicity, and AlphaMissense scores. Local tooling (Ensembl VEP CLI with a 50 GB cache, ANNOVAR with annotation databases) is reproducible but heavyweight; public REST endpoints (Ensembl VEP, OncoKB, ClinVar, gnomAD) are programmatically accessible and free, with rate limits well above what a 36-sample or 1000-sample cohort needs.

## Decision

Phase 1 annotation goes through a pluggable `Annotator` interface with multiple backends selected by a CLI flag. Defaults favor API-backed implementations over local-cache implementations.

The same pattern applies to future swappable services in other phases (structure predictors in phase 2, dependency sources in phase 3, drug catalogs in phase 4). Each phase that grows multiple backends gets its own sibling package: `bin/annotators/`, eventually `bin/structure_predictors/`, etc.

## Consequences

- The first concrete realization: `bin/annotators/{base, curated_vcf, vep_rest}.py` with `--annotator {curated, vep_rest}` on `bin/phase1_annotate.py`.
- Defaults are explicit per-script. Toy fixtures (used by the scorecard) default to `curated` so the case-study scorecard remains fast and offline. Real-cohort runs (`bin/run_target_rt.py`) explicitly pass `--annotator vep_rest`.
- API-backed backends cache to `data/<service>_cache/` (already gitignored). Cache invalidation is `rm -rf` of the directory. No reference data is committed to the repo.
- CI does not hit the network. Tests mock the HTTP layer.
- Local-CLI backends (e.g., a hypothetical `bin/annotators/vep_local.py` that shells out to a containerized VEP) are explicitly allowed but only added when an API-backed backend cannot satisfy a real requirement (rate-limit pain, plugin needs, regulatory reproducibility).

## Addendum (v0.13): chained annotators

The `--annotator` flag accepts a comma-separated list (e.g., `vep_rest,oncokb`). The first member is a primary backend that produces a complete `VariantAnnotation`; subsequent members may either:

- override the prior annotation (calling `annotate_batch` again), or
- augment it via the `AugmentingAnnotator.annotate_with_prior` protocol, which receives both the original variants and the prior pass's annotations.

OncoKB is the first augmenting backend. It consumes `gene` + `hgvsp_short` from the prior pass and adds `oncogenic` + `mutation_effect` fields to the annotation. `bin/phase1_annotate.py:classify_snv` honors these new fields by lifting any protein-altering missense to DRIVER when OncoKB calls it Oncogenic / Likely Oncogenic, even if our hand-curated hotspot list does not include the residue.

Authentication and disease parameterization for augmenting backends:

- OncoKB token comes from `ONCOKB_TOKEN` env var (or constructor arg for tests). When the token is missing, the OncoKB pass logs a warning once and passes prior annotations through unchanged. Real-cohort runs continue to produce VEP-only output.
- `--oncokb-tumor-type` accepts an OncoTree code or generic name. Default `RMS`; subtype-specific `ERMS` / `ARMS` / `PLRMS` available for callers that want narrower context.

Environment-driven enablement: `bin/run_target_rt.py` auto-appends `,oncokb` to the chain when `ONCOKB_TOKEN` is set, so cohort runs upgrade automatically once a user obtains an academic token. No code edit required to opt in.

## Alternatives considered

- **Single-backend phase 1, no abstraction.** Rejected: the very next thing the pipeline needs (OncoKB layer in v0.13) would force the same refactor with worse motivation than building it now.
- **Local VEP CLI as the default.** Rejected: 50 GB cache; CI can't run it; iteration speed suffers; doesn't generalize to OncoKB/ClinVar/AlphaMissense, which are API-first by design.
- **All annotation upstream, in fetcher scripts.** Rejected: couples annotation choice to data-source choice. The toy fixtures don't have a fetcher; future cohorts (MCI) will have a totally different fetcher; the annotator should be a phase-1 concern, not a fetch concern.
- **Standalone OncoKB phase, not chained with VEP.** Rejected for v0.13: phase 1 is the natural home for any per-variant annotation, and the chain pattern generalizes to ClinVar / AlphaMissense / future backends without inventing new pipeline phases.
