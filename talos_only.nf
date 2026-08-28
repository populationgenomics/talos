#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { TALOS } from './nextflow/talos'
include { SV_ANNOTATION } from './nextflow/sv_annotation'

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
	ch_gff = channel.fromPath(params.ensembl_gff, checkIfExists: true).first()
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
			files("${workflow.outputDir.toUriString()}/${row.cohort}_outputs/*.mt", type: 'dir'),
			file(row.pedigree, checkIfExists: true),
			file(row.config, checkIfExists: true),
			// optional columns - an empty/absent cell becomes [], which stages nothing and is falsy in every downstream truthiness check
			row.history ? file(row.history, checkIfExists: true) : [],
			row.ext_ids ? file(row.ext_ids, checkIfExists: true) : [],
			row.seqr_map ? file(row.seqr_map, checkIfExists: true) : [],
			row.mito ? file(row.mito, checkIfExists: true) : [],
		) }

	// the SV path is entirely optional, and only wired up if the input TSV declares an `sv` column.
	// this is checked eagerly rather than per-row so that a cohort with no SV data never requires the
	// SV reference files to be present at all
	def sv_requested = file(params.input_tsv).withReader { handle -> handle.readLine() }.tokenize('\t').contains('sv')

	ch_sv_annotated = channel.empty()
	ch_sv_fresh = channel.empty()

	if (sv_requested) {
		ch_sv_inputs = channel.fromPath(params.input_tsv)
			.splitCsv(header: true, sep: '\t')
			.map { row -> tuple(
				row.cohort,
				row.sv ? file(row.sv, checkIfExists: true) : [],
				file(row.config, checkIfExists: true),
			) }

		SV_ANNOTATION(
			ch_ref_genome,
			ch_sv_inputs,
		)

		ch_sv_annotated = SV_ANNOTATION.out.annotated
		ch_sv_fresh = SV_ANNOTATION.out.fresh
	}

	TALOS(
		ch_mane,
		ch_gff,
        ch_ref_genome,
		ch_inputs,
		ch_sv_annotated,
	)

	publish:
    	html = TALOS.out.html
		json = TALOS.out.json
		panelapp = TALOS.out.panelapp
		labelled_sv = TALOS.out.labelled_sv
		sv_annotated = ch_sv_fresh
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
	labelled_sv {
		path { id, _labelled_sv, _labelled_sv_idx -> "${id}_outputs" }
	}
	sv_annotated {
		path { id, _vcf, _vcf_idx -> "${id}_outputs" }
	}
}
