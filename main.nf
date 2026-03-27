#!/usr/bin/env nextflow

/*
    Talos Unified Workflow
    ======================

    This is the main entry point for the Talos pipeline. It orchestrates two distinct workflows:

    1. ANNOTATION: Annotates VCF files with the prepared data.
    2. TALOS: Runs the core Talos analysis/filtering/reporting.

    Usage:
    nextflow run nextflow/main.nf --input_tsv [path] [other params...]
*/

nextflow.enable.dsl=2

// Import specific workflows
include { ANNOTATION } from './nextflow/annotation'
include { TALOS } from './nextflow/talos'


workflow {
	main:
	if (file(workflow.outputDir).simpleName == file(params.processed_annotations).simpleName) {
    	println "Output Directory (${workflow.outputDir}) is probably not set correctly, use config or `-output-dir`"
		exit 1
    }

	if (!file(params.mane_json).exists()) {
		println "MANE JSON not available, please run the Talos Prep workflow (--entry preparation)"
		exit 1
	}

	if (!params.input_tsv) {
		println "Required --input_tsv argument not provided"
		exit 1
	}

	ch_gff = Channel.fromPath(params.ensembl_gff, checkIfExists: true)
	ch_ref_genome = Channel.fromPath(params.ref_genome, checkIfExists: true)
	ch_mane = Channel.fromPath(params.mane_json, checkIfExists: true)

	ch_inputs = Channel.fromPath(params.input_tsv)
		.splitCsv(header: true, sep: '\t')
		.map { row -> tuple(row.cohort, row.path, row.type) }

	ch_talos_inputs = Channel.fromPath(params.input_tsv)
		.splitCsv(header: true, sep: '\t')
		.map { row -> tuple(
			row.cohort,
			row.path,
			row.type,
			file(row.pedigree, checkIfExists: true),
			file(row.config, checkIfExists: true),
			file(row.history, checkIfExists: true),
			file(row.ext_ids, checkIfExists: true),
			file(row.seqr_map, checkIfExists: true),
		) }

	ANNOTATION(
		ch_gff,
		ch_mane,
		ch_ref_genome,
		ch_inputs,
	)

	ch_talos_combined = ANNOTATION.out.mts
		.join(ch_talos_inputs)
		.map { cohort, mts, inpath, intype, pedigree, config, history, ext, seqr -> tuple(cohort, mts, pedigree, config, history, ext, seqr) }

	TALOS(
		ch_mane,
		ch_talos_combined,
	)

	publish:
		mts = ANNOTATION.out.mts
    	html = TALOS.out.html
		json = TALOS.out.json
		panelapp = TALOS.out.panelapp
}

output {
	mts {
		path { id, mts -> "${id}_outputs" }
	}
	html {
		path { id, html -> "${id}_outputs" }
	}
	json {
		path { id, json -> "${id}_outputs" }
	}
	panelapp {
		path { id, panelapp -> "${id}_outputs" }
	}
}
