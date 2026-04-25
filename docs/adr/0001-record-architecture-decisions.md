# ADR 0001: Record Architecture Decisions

Date: 2026-04-25

## Status

Accepted.

## Context

This pipeline involves many non-obvious choices: which orchestrator (Nextflow vs Snakemake), which structure predictor (AlphaFold DB vs Boltz vs Chai), which drug database hierarchy (OncoKB vs MTP vs DGIdb first), how to weight confidence components, when to spend GPU time. Without a record, those choices become invisible institutional knowledge.

## Decision

We use lightweight Architecture Decision Records (ADRs), one per non-trivial choice, stored under `docs/adr/`. Each ADR is short (under one page), dated, and immutable once accepted. Subsequent reversals get a new ADR that supersedes the old one.

Format:

```
# ADR NNNN: Title
Date: YYYY-MM-DD
Status: Proposed | Accepted | Superseded by ADR XXXX
Context: ...
Decision: ...
Consequences: ...
```

## Consequences

Future contributors (including future-self) can audit the rationale behind each design choice without reading every commit message. Reviewers can challenge a decision by writing the counter-ADR, not by relitigating it conversationally.
