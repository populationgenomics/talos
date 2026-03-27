#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { TALOS } from './nextflow/talos'

workflow {
	main:
    if (file(workflow.outputDir).simpleName == file(params.processed_annotations).simpleName) {
    	println "Output Directory (${workflow.outputDir}) is probably not set correctly, use config or `-output-dir`"
		exit 1
    }

	if (!params.input_tsv) {
		println "Required --input_tsv argument not provided"
		exit 1
	}
	ch_mane = Channel.fromPath(params.mane_json, checkIfExists: true)
	ch_inputs = Channel.fromPath(params.input_tsv)
		.splitCsv(header: true, sep: '\t')
		.map { row -> tuple(
			row.cohort,
			Channel.fromPath("${outputDir}/${row.cohort}_outputs/*.mt", checkIfExists: true),
			file(row.pedigree, checkIfExists: true),
			file(row.config, checkIfExists: true),
			file(row.history, checkIfExists: true),
			file(row.ext_ids, checkIfExists: true),
			file(row.seqr_map, checkIfExists: true),
		) }

	TALOS(
		ch_mane,
		ch_inputs,
	)

	publish:
    	html = TALOS.out.html
		json = TALOS.out.json
		panelapp = TALOS.out.panelapp
}

output {
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
