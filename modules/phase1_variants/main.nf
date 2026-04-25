// modules/phase1_variants/main.nf
// Phase 1: Variant discovery and annotation.
// Skeleton only. Implementation begins at Step 2.1 of PIPELINE_STANDUP_PLAN.md.

nextflow.enable.dsl=2

process PHASE1_ANNOTATE {
    tag "phase1:${meta.sample_id}"
    container "ghcr.io/mzobeck/rms-isp-phase1:0.1.0"
    publishDir "${params.outdir}/phase1", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.sample_id}.annotated.tsv"), emit: annotated
    tuple val(meta), path("${meta.sample_id}.annotation_log.txt"), emit: log

    script:
    """
    echo "PHASE1_ANNOTATE not yet implemented." > ${meta.sample_id}.annotation_log.txt
    echo "See PIPELINE_STANDUP_PLAN Phase 2 (Steps 2.1-2.4) for implementation plan." >> ${meta.sample_id}.annotation_log.txt
    head -1 ${vcf} > ${meta.sample_id}.annotated.tsv || true
    """
}
