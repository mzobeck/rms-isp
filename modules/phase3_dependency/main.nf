// modules/phase3_dependency/main.nf
// Phase 3: Expression and dependency integration.
// v0.1: joins against assets/depmap_rms_summary.tsv (curated PLACEHOLDER).
// v0.2 will swap to a real DepMap 24Q2+ pull and add OpenPedCan expression z-scores.

nextflow.enable.dsl=2

process PHASE3_DEPENDENCY {
    tag "phase3:${sample_id}"
    publishDir "${params.outdir}/phase3", mode: 'copy'

    input:
    tuple val(sample_id), path(structured)
    path depmap_summary
    path expression_summary
    val subtype

    output:
    tuple val(sample_id), path("${sample_id}.phase3.tsv"), emit: dependency

    script:
    def expr_arg = expression_summary.name == 'NO_EXPR' ? '' : "--expression ${expression_summary}"
    """
    phase3_dependency.py \\
        --in ${structured} \\
        --depmap ${depmap_summary} \\
        ${expr_arg} \\
        --subtype ${subtype} \\
        --out ${sample_id}.phase3.tsv
    """
}
