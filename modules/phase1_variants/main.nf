// modules/phase1_variants/main.nf
// Phase 1: Variant discovery and annotation.
// v0.1: joins VCF against assets/targets_kb.tsv to call DRIVER / VUS / PASSENGER / OFF_TARGET.
// Real implementation will run VEP + OncoKB; the output column contract is stable.

nextflow.enable.dsl=2

process PHASE1_ANNOTATE {
    tag "phase1:${sample_id}"
    publishDir "${params.outdir}/phase1", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf)
    path targets_kb

    output:
    tuple val(sample_id), path("${sample_id}.phase1.tsv"), emit: annotated

    script:
    """
    phase1_annotate.py \\
        --vcf ${vcf} \\
        --targets-kb ${targets_kb} \\
        --sample-id ${sample_id} \\
        --out ${sample_id}.phase1.tsv
    """
}
