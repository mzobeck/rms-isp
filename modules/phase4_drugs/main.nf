// modules/phase4_drugs/main.nf
// Phase 4: Drug matching.
// v0.4: unions assets/drug_target_map.tsv (curated, primary; mechanism strings
// drive the case-study scorecard) with assets/dgidb_drugs.tsv (auto-generated
// from DGIdb GraphQL; refreshed by bin/fetch_dgidb.py). Curated wins on
// (gene, drug) collisions.

nextflow.enable.dsl=2

process PHASE4_DRUGS {
    tag "phase4:${sample_id}"
    publishDir "${params.outdir}/phase4", mode: 'copy'

    input:
    tuple val(sample_id), path(dependency)
    path drug_map
    path drug_map_extra

    output:
    tuple val(sample_id), path("${sample_id}.phase4.tsv"), emit: drugs

    script:
    def extra_arg = drug_map_extra.name == 'NO_EXTRA' ? '' : "--drug-map-extra ${drug_map_extra}"
    """
    phase4_drugs.py \\
        --in ${dependency} \\
        --drug-map ${drug_map} \\
        ${extra_arg} \\
        --out ${sample_id}.phase4.tsv
    """
}
