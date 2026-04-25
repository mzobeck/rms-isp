// modules/phase2_structure/main.nf
// Phase 2: Structural reference attachment.
// v0.1: AlphaFold DB URL + heuristic structural-confidence score.
// v0.2 will add Boltz-1 / Chai-1 mutant prediction and AutoDock Vina docking on a shortlist.

nextflow.enable.dsl=2

process PHASE2_STRUCTURE {
    tag "phase2:${sample_id}"
    publishDir "${params.outdir}/phase2", mode: 'copy'

    input:
    tuple val(sample_id), path(annotated)

    output:
    tuple val(sample_id), path("${sample_id}.phase2.tsv"), emit: structured

    script:
    """
    phase2_structure.py \\
        --in ${annotated} \\
        --out ${sample_id}.phase2.tsv
    """
}
