// modules/phase4_drugs/main.nf
// Phase 4: Drug matching.
// v0.1: joins against assets/drug_target_map.tsv (curated subset).
// v0.2 will add live DGIdb v5, OpenTargets, and ClinicalTrials.gov queries.

nextflow.enable.dsl=2

process PHASE4_DRUGS {
    tag "phase4:${sample_id}"
    publishDir "${params.outdir}/phase4", mode: 'copy'

    input:
    tuple val(sample_id), path(dependency)
    path drug_map

    output:
    tuple val(sample_id), path("${sample_id}.phase4.tsv"), emit: drugs

    script:
    """
    phase4_drugs.py \\
        --in ${dependency} \\
        --drug-map ${drug_map} \\
        --out ${sample_id}.phase4.tsv
    """
}
