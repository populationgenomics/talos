#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
    Talos Unified Workflow
    ======================

    This is the main entry point for the Talos pipeline. It orchestrates two distinct workflows:

    1. ANNOTATION: Annotates VCF file(s) with the prepared data.
    2. TALOS: Runs the core Talos analysis/filtering/reporting.

    Usage:
    nextflow run nextflow/main.nf --input_tsv [path] [other params...]
*/

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

	ch_gff = channel.fromPath(params.ensembl_gff, checkIfExists: true).first()
	ch_ref_genome = channel.fromPath(params.ref_genome, checkIfExists: true).first()
	ch_mane = channel.fromPath(params.mane_json, checkIfExists: true).first()

	ch_inputs = channel.fromPath(params.input_tsv)
		.splitCsv(header: true, sep: '\t')
		.map { row -> tuple(row.cohort, row.path, row.type) }

	ch_talos_inputs = channel.fromPath(params.input_tsv)
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
			file(row.mito, checkIfExists: true),
		) }

	ANNOTATION(
		ch_gff,
		ch_mane,
		ch_ref_genome,
		ch_inputs,
	)

	ch_talos_combined = ANNOTATION.out.mts
		.join(ch_talos_inputs)
		.map { cohort, mts, _inpath, _intype, pedigree, config, history, ext, seqr, mito -> tuple(cohort, mts, pedigree, config, history, ext, seqr, mito) }

	TALOS(
		ch_mane,
		ch_gff,
		ch_ref_genome,
		ch_talos_combined,
	)

	publish:
		mts = ANNOTATION.out.mts
    	html = TALOS.out.html
		json = TALOS.out.json
		labelled = TALOS.out.labelled
		panelapp = TALOS.out.panelapp
}

output {
	mts {
		path { id, _mts -> "${id}_outputs" }
	}
	html {
		path { id, _html -> "${id}_outputs" }
	}
	json {
		path { id, _json -> "${id}_outputs" }
	}
	panelapp {
		path { id, _panelapp -> "${id}_outputs" }
	}
	labelled {
		path { id, _labelled, _labelled_idx -> "${id}_outputs" }
	}
}
