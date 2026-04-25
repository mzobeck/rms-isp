#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * RMS-ISP: Rhabdomyosarcoma In Silico Pipeline
 * v0.1.0-pilot: end-to-end variant -> annotated -> structure -> dependency -> drugs -> ranked report.
 * See plans/pilot_ccdi_mci/rms_translational_pilot_project.md for the build plan.
 */

include { PHASE1_ANNOTATE  } from './modules/phase1_variants/main.nf'
include { PHASE2_STRUCTURE } from './modules/phase2_structure/main.nf'
include { PHASE3_DEPENDENCY } from './modules/phase3_dependency/main.nf'
include { PHASE4_DRUGS     } from './modules/phase4_drugs/main.nf'
include { PHASE5_SCORE     } from './modules/phase5_scoring/main.nf'

// All param defaults live in nextflow.config; CLI flags override.
params.pipeline_version = workflow.manifest.version

log.info """
=========================================================
RMS-ISP  v${params.pipeline_version}
---------------------------------------------------------
input          : ${params.input}
sample_id      : ${params.sample_id}
subtype        : ${params.subtype}
outdir         : ${params.outdir}
targets_kb     : ${params.targets_kb}
depmap_summary : ${params.depmap_summary}
drug_map       : ${params.drug_map}
profile        : ${workflow.profile}
=========================================================
""".stripIndent()

workflow {
    vcf_ch         = Channel.fromPath(params.input, checkIfExists: true)
    targets_kb_ch  = Channel.value(file(params.targets_kb,     checkIfExists: true))
    depmap_ch      = Channel.value(file(params.depmap_summary, checkIfExists: true))
    drug_map_ch    = Channel.value(file(params.drug_map,       checkIfExists: true))

    sample_vcf_ch  = vcf_ch.map { v -> tuple(params.sample_id, v) }

    PHASE1_ANNOTATE(sample_vcf_ch, targets_kb_ch)
    PHASE2_STRUCTURE(PHASE1_ANNOTATE.out.annotated)
    PHASE3_DEPENDENCY(PHASE2_STRUCTURE.out.structured, depmap_ch, params.subtype)
    PHASE4_DRUGS(PHASE3_DEPENDENCY.out.dependency, drug_map_ch)

    p5_input_ch = PHASE4_DRUGS.out.drugs
        .combine(vcf_ch)
        .map { sid, drugs, vcf -> tuple(sid, drugs, vcf) }

    PHASE5_SCORE(p5_input_ch, params.subtype, params.pipeline_version)

    PHASE5_SCORE.out.report.view { sid, md ->
        "report ready: ${sid} -> ${md}"
    }
}
