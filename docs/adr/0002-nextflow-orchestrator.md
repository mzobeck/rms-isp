# ADR 0002: Nextflow as Pipeline Orchestrator

Date: 2026-04-25

## Status

Accepted.

## Context

Two viable choices: Nextflow and Snakemake. Both are mature, both support Docker/Singularity, both run on Slurm and cloud batch.

## Decision

Nextflow.

Reasoning:

- Native, first-class support for AWS Batch, GCP Batch, and Slurm without rewrite.
- nf-core community provides production-grade module patterns directly applicable to genomic workflows.
- PedcBioPortal, OpenPedCan, and many CCDI-affiliated pipelines use Nextflow; alignment lowers the cost of cross-pipeline reuse.
- Seqera Platform (Tower) provides drop-in provenance and monitoring if we choose to adopt it later.

Snakemake's primary advantages (simpler Python-native syntax, smaller learning curve) do not outweigh the cloud and ecosystem alignment.

## Consequences

We commit to DSL2 modules. Process-level idioms follow nf-core conventions where applicable. Snakemake interoperability is not a goal.
