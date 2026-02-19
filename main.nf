#!/usr/bin/env nextflow

/*
    Talos Unified Workflow
    ======================

    This is the main entry point for the Talos pipeline. It orchestrates two distinct workflows:

    1. ANNOTATION: Annotates VCF files with the prepared data.
    2. TALOS: Runs the core Talos analysis/filtering/reporting.

    Usage:
    nextflow run nextflow/main.nf -entry [ANNOTATION, TALOS] [other params...]
*/

// Import specific workflows
include { ANNOTATION } from './nextflow/annotation'
include { TALOS } from './nextflow/talos'

workflow TALOS_ONLY {
	ch_mane = Channel.fromPath(params.mane_json, checkIfExists: true)
	ch_mts = Channel.fromPath("${params.cohort_output_dir}/*.mt", type: 'dir', checkIfExists: true).collect()
	TALOS(
		ch_mane,
		ch_mts,
	)
}

workflow {
	ch_gff = Channel.fromPath(params.ensembl_gff, checkIfExists: true)
	ch_ref_genome = Channel.fromPath(params.ref_genome, checkIfExists: true)

	if (!file(params.mane_json).exists()) {
		println "MANE JSON not available, please run the Talos Prep workflow (--entry preparation)"
		exit 1
	}

	ch_mane = Channel.fromPath(params.mane_json, checkIfExists: true)

	ANNOTATION(
		ch_gff,
		ch_mane,
		ch_ref_genome,
	)

	TALOS(
		ch_mane,
		ANNOTATION.out.mts,
	)
}
