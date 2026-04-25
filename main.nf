#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * RMS-ISP: Rhabdomyosarcoma In Silico Pipeline
 * Top-level workflow stub. Phases will be wired in as modules complete.
 */

params.input          = params.input          ?: "${projectDir}/tests/data/toy_patient.vcf"
params.outdir         = params.outdir         ?: "${projectDir}/results"
params.targets_yaml   = params.targets_yaml   ?: "${projectDir}/conf/targets.yaml"

log.info """
=========================================================
RMS-ISP  v0.0.1  (scaffold)
---------------------------------------------------------
input         : ${params.input}
outdir        : ${params.outdir}
targets_yaml  : ${params.targets_yaml}
profile       : ${workflow.profile}
=========================================================
""".stripIndent()

process HELLO {
    tag "hello"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path vcf

    output:
    path "hello.txt"

    script:
    """
    echo "RMS-ISP scaffold OK." > hello.txt
    echo "Received input: ${vcf.name}" >> hello.txt
    echo "Pipeline scaffold smoke test passed." >> hello.txt
    """
}

workflow {
    Channel.fromPath(params.input, checkIfExists: true) | HELLO
}
