// modules/phase5_scoring/main.nf
// Phase 5: Confidence scoring + patient-level markdown report.
// v0.1: weighted-sum confidence per pilot plan §4 Aim 2 Phase 5; expression
// component is deferred and contributes 0.

nextflow.enable.dsl=2

process PHASE5_SCORE {
    tag "phase5:${sample_id}"
    publishDir "${params.outdir}/phase5", mode: 'copy'

    input:
    tuple val(sample_id), path(drugs), path(vcf)
    val subtype
    val pipeline_version

    output:
    tuple val(sample_id), path("${sample_id}.phase5.tsv"), emit: scored
    tuple val(sample_id), path("${sample_id}.report.md"),  emit: report

    script:
    """
    phase5_score.py \\
        --in ${drugs} \\
        --vcf ${vcf} \\
        --sample-id ${sample_id} \\
        --subtype ${subtype} \\
        --pipeline-version ${pipeline_version} \\
        --out-tsv ${sample_id}.phase5.tsv \\
        --out-md  ${sample_id}.report.md
    """
}
