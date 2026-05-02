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
	ch_mane = channel.fromPath(params.mane_json, checkIfExists: true).first()
	ch_ref_genome = channel.fromPath(params.ref_genome, checkIfExists: true).first()

	/*
	NextFlow doesn't like the use of files("...", type: 'dir') here, as it scans for files at DAG setup time instead of
	each time a new Stage is reached. This is intentional, as these files are not created by this workflow
	*/
	ch_inputs = channel.fromPath(params.input_tsv)
		.splitCsv(header: true, sep: '\t')
		.map { row -> tuple(
			row.cohort,
			files("${workflow.outputDir}/${row.cohort}_outputs/*.mt", type: 'dir'),
			file(row.pedigree, checkIfExists: true),
			file(row.config, checkIfExists: true),
			file(row.history, checkIfExists: true),
			file(row.ext_ids, checkIfExists: true),
			file(row.seqr_map, checkIfExists: true),
		) }

	TALOS(
		ch_mane,
		ch_gff,
        ch_ref_genome,
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
