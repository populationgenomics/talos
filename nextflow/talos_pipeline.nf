#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
The shared Talos pipeline: ANNOTATION -> SV_ANNOTATION -> TALOS, for any number of cohorts.

Entrypoints differ only in how they build the cohort rows and decide where each cohort's outputs go (see
docs/adr/0001). Everything from here down is layout-agnostic: each row carries its own output directories,
and nothing below this workflow builds an output path.
*/

include { ANNOTATION } from './annotation'
include { SV_ANNOTATION } from './sv_annotation'
include { TALOS } from './talos'

// an optional input column: the checked file, or [] when the cell is empty/absent
def optFile(path) {
	path ? file(path, checkIfExists: true) : []
}

workflow TALOS_PIPELINE {
	take:
		// one Map per cohort: every input TSV column (cohort, path, type, pedigree, config, and the optional
		// history, ext_ids, seqr_map, mito, sv, sv_annotated, str), plus annotated_dir and analysis_dir - the
		// String URIs this cohort's annotation products and analysis results are published to - and optionally
		// report_dir, where the HTML report is published (default: analysis_dir)
		ch_cohorts_raw

	main:
	if (!file(params.mane_json).exists()) {
		println "MANE JSON not available, please run the Talos Prep workflow (--entry preparation)"
		exit 1
	}

	if (!params.ensembl_gff) {
		println "params.ensembl_gff (${params.ensembl_gff}) is not available, please re-run the input file download script & preparation.nf workflow"
		exit 1
	}

	ch_gff = channel.fromPath(params.ensembl_gff, checkIfExists: true).first()
	ch_ref_genome = channel.fromPath(params.ref_genome, checkIfExists: true).first()
	ch_mane = channel.fromPath(params.mane_json, checkIfExists: true).first()

	// `sv` is a raw joint-call to annotate, `sv_annotated` an already annotated one (e.g. GATK-SV output) that
	// skips SV_ANNOTATION. Both on one row is ambiguous - which one should be analysed?
	// Every file column is resolved (and checked) once here; the channels below are plain projections of this row.
	// Optional columns become [] when absent, which stages nothing and is falsy in every downstream truthiness check
	ch_cohorts = ch_cohorts_raw.map { row ->
		if (row.sv && row.sv_annotated) {
			error("Cohort ${row.cohort} sets both `sv` (${row.sv}) and `sv_annotated` (${row.sv_annotated}), supply one or the other")
		}
		row + [
			pedigree: file(row.pedigree, checkIfExists: true),
			config: file(row.config, checkIfExists: true),
			history: optFile(row.history),
			ext_ids: optFile(row.ext_ids),
			seqr_map: optFile(row.seqr_map),
			mito: optFile(row.mito),
			sv: optFile(row.sv),
			// STR VCFs (from stripy_json_to_vcf) need no annotation, they go straight to ValidateMOI
			str: optFile(row.str),
			report_dir: row.report_dir ?: row.analysis_dir,
			// whether TALOS should wait for a filtered SV VCF for this cohort
			has_sv: (row.sv || row.sv_annotated) as boolean,
		]
	}

	ch_inputs = ch_cohorts.map { row -> tuple(row.cohort, row.path, row.type, row.pedigree, row.annotated_dir) }

	// per-cohort metadata, one row per cohort as [cohort, row] - the annotated shards are carried separately
	ch_talos_meta = ch_cohorts.map { row -> tuple(row.cohort, row) }

	// SV annotation is wired up for every cohort, and decided per row: a cohort with no `sv` VCF carries [],
	// and drops out in SV_ANNOTATION's no_sv branch. The SV reference data is therefore always required
	ch_sv_inputs = ch_cohorts.map { row -> tuple(row.cohort, row.sv, row.config, row.pedigree, row.annotated_dir) }

	SV_ANNOTATION(
		ch_ref_genome,
		ch_sv_inputs,
	)

	// pre-annotated SV VCFs bypass SV_ANNOTATION, joining its output in the same [cohort, vcf, tbi] shape.
	// Their AF fields aren't renamed, so the cohort config's RunSvFiltering.gnomad_population must match the
	// VCF's own field prefix
	ch_sv_preannotated = ch_cohorts
		.filter { row -> row.sv_annotated }
		.map { row -> tuple(
			row.cohort,
			file(row.sv_annotated, checkIfExists: true),
			file("${row.sv_annotated}.tbi", checkIfExists: true),
		) }

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
		ANNOTATION.out.shards,
		SV_ANNOTATION.out.annotated.mix(ch_sv_preannotated),
	)

	emit:
		// one entry per cohort: [row, published_uri] - the cohort row Map (as resolved above), and the String URI
		// its results JSON / HTML report is published to. The published copy doesn't exist until publishDir has
		// finished copying it, so a consumer that reads it must wait for it
		json = TALOS.out.json
			.join(ch_talos_meta)
			.map { _cohort, results_json, row -> tuple(row, "${row.analysis_dir}/${results_json.name}".toString()) }
		html = TALOS.out.html
			.join(ch_talos_meta)
			.map { _cohort, report_html, row -> tuple(row, "${row.report_dir}/${report_html.name}".toString()) }
}
