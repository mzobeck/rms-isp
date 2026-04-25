// modules/phase1_variants/main.nf
// Phase 1: Variant + CNA + fusion annotation.
// v0.2: unified event table (snv | cna | fusion) joined against assets/targets_kb.tsv.
// Real implementation will run VEP + OncoKB for SNVs, GISTIC2/cnvkit/FACETS for CNAs,
// and STAR-Fusion / Arriba for fusions; the output column contract is stable.

nextflow.enable.dsl=2

process PHASE1_ANNOTATE {
    tag "phase1:${sample_id}"
    publishDir "${params.outdir}/phase1", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(cna), path(fusion)
    path targets_kb

    output:
    tuple val(sample_id), path("${sample_id}.phase1.tsv"), emit: annotated

    script:
    """
    phase1_annotate.py \\
        --vcf ${vcf} \\
        --cna ${cna} \\
        --fusion ${fusion} \\
        --targets-kb ${targets_kb} \\
        --sample-id ${sample_id} \\
        --out ${sample_id}.phase1.tsv
    """
}
