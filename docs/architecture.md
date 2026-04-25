# Architecture

Placeholder. To be written at Phase 1 close.

## Intended contents

- Overall data flow (variant input -> annotated -> structural -> dependency -> drug -> scored report).
- Module boundaries: each phase is independently runnable, with defined input and output schemas.
- Provenance and reproducibility model: per-run manifests recording pipeline version, container digests, parameters.
- Configuration model: profiles for laptop, slurm, aws, gpu.
- Testing strategy: unit tests per module, integration test on toy patient, golden output diffs.

See `../proj_management/PIPELINE_STANDUP_PLAN.md` for the build plan.
