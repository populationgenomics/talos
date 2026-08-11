#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
    Talos Unified Workflow
    ======================

    This is the main entry point for the Talos pipeline. It orchestrates two distinct workflows:

    1. ANNOTATION: Annotates VCF file(s) with the prepared data.
    2. SV_ANNOTATION: Annotates a joint-called SV VCF, if the input TSV carries an `sv` column.
    3. TALOS: Runs the core Talos analysis/filtering/reporting.

    Usage:
    nextflow run nextflow/main.nf --input_tsv [path] [other params...]
*/

// Import specific workflows
include { ANNOTATION } from './nextflow/annotation'
include { SV_ANNOTATION } from './nextflow/sv_annotation'
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

	if (!params.ensembl_gff) {
		println "params.ensembl_gff (${params.ensembl_gff}) is not available, please re-run the input file download script & preparation.nf workflow"
		exit 1
	}

	ch_gff = channel.fromPath(params.ensembl_gff, checkIfExists: true).first()
	ch_ref_genome = channel.fromPath(params.ref_genome, checkIfExists: true).first()
	ch_mane = channel.fromPath(params.mane_json, checkIfExists: true).first()

	ch_inputs = channel.fromPath(params.input_tsv)
		.splitCsv(header: true, sep: '\t')
		.map { row -> tuple(row.cohort, row.path, row.type) }

	// per-cohort metadata, one row per cohort - the annotated shards are carried separately
	ch_talos_meta = channel.fromPath(params.input_tsv)
		.splitCsv(header: true, sep: '\t')
		.map { row -> tuple(
			row.cohort,
			file(row.pedigree, checkIfExists: true),
			file(row.config, checkIfExists: true),
			file(row.history ?: "${projectDir}/nextflow/assets/NO_HISTORY", checkIfExists: true),
			file(row.ext_ids ?: "${projectDir}/nextflow/assets/NO_FILE", checkIfExists: true),
			file(row.seqr_map ?: "${projectDir}/nextflow/assets/NO_SEQR_FILE", checkIfExists: true),
			file(row.mito ?: "${projectDir}/nextflow/assets/NO_MITO", checkIfExists: true),
		) }

	// the SV path is entirely optional, and only wired up if the input TSV declares an `sv` column.
	// this is checked eagerly rather than per-row so that a cohort with no SV data never requires the
	// SV reference files to be present at all
	def sv_requested = file(params.input_tsv).withReader { handle -> handle.readLine() }.tokenize('\t').contains('sv')

	ch_sv_annotated = channel.empty()

	if (sv_requested) {
		ch_sv_inputs = channel.fromPath(params.input_tsv)
			.splitCsv(header: true, sep: '\t')
			.map { row -> tuple(
				row.cohort,
				file(row.sv ?: "${projectDir}/nextflow/assets/NO_SV", checkIfExists: true),
				file(row.config, checkIfExists: true),
			) }

		SV_ANNOTATION(
			ch_ref_genome,
			ch_sv_inputs,
		)

		ch_sv_annotated = SV_ANNOTATION.out.annotated
	}

	ANNOTATION(
		ch_gff,
		ch_ref_genome,
		ch_inputs,
	)

	TALOS(
		ch_mane,
		ch_gff,
		ch_ref_genome,
		ch_talos_meta,
		ANNOTATION.out.vcfs,
		ch_sv_annotated,
	)

	publish:
		annotated = ANNOTATION.out.vcfs
    	html = TALOS.out.html
		json = TALOS.out.json
		labelled = TALOS.out.labelled
		labelled_sv = TALOS.out.labelled_sv
		panelapp = TALOS.out.panelapp
		sv_annotated = ch_sv_annotated
}

output {
	annotated {
		path { id, _vcf -> "${id}_outputs" }
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
	labelled_sv {
		path { id, _labelled_sv, _labelled_sv_idx -> "${id}_outputs" }
	}
	sv_annotated {
		path { id, _vcf, _vcf_idx -> "${id}_outputs" }
	}
}
