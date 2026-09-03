#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { TALOS } from './nextflow/talos'
include { SV_ANNOTATION } from './nextflow/sv_annotation'

// analysis outputs are published to a per-run dated directory, ${cohort}_analysis_YYYYMMDD -
// annotation products go to the undated ${cohort}_annotated, which is reused across cycles
def runDate() {
    workflow.start.format(java.time.format.DateTimeFormatter.ofPattern('yyyyMMdd'))
}

// build the content of a proposed input TSV for the next reanalysis cycle: this run's input TSV
// verbatim, except each cohort's `history` cell now points at the results JSON this run produced.
// resultsByCohort maps cohort -> published results JSON path; cohorts absent from it keep whatever
// history they were given
def nextInputTsv(input_tsv, resultsByCohort) {
    def lines = input_tsv.readLines().findAll { line -> line.trim() }
    // split with a negative limit - trailing empty cells matter here, they are the optional columns
    def header = lines[0].split('\t', -1) as List
    if (!header.contains('history')) {
        header << 'history'
    }
    def cohort_idx = header.indexOf('cohort')
    def history_idx = header.indexOf('history')

    def rows = [header.join('\t')]
    lines.drop(1).each { line ->
        def raw = line.split('\t', -1) as List
        // rows commonly omit trailing empty cells altogether, so pad out to the header width
        def cells = raw.size() < header.size() ? raw + [''] * (header.size() - raw.size()) : raw
        def cohort = cells[cohort_idx]
        if (resultsByCohort.containsKey(cohort)) {
            cells[history_idx] = resultsByCohort[cohort]
        }
        rows << cells.join('\t')
    }
    return rows.join('\n') + '\n'
}

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

	// a ready-to-use input TSV for the next reanalysis cycle - identical to this run's input, with
	// `history` repointed at the results JSON published for that cohort by this run. The paths are
	// the published destinations, which don't exist until the run completes, so this is a proposal
	// to be reviewed, not an output consumed by anything in this run
	ch_next_input_tsv = TALOS.out.json
		.map { cohort, results_json -> tuple(cohort, "${workflow.outputDir}/${cohort}_analysis_${runDate()}/${results_json.name}".toString()) }
		.collect(flat: false)
		.map { pairs -> nextInputTsv(file(params.input_tsv), pairs.collectEntries { pair -> pair }) }
		.collectFile(name: "talos_input_${runDate()}.tsv")

	publish:
		next_input_tsv = ch_next_input_tsv
    	html = TALOS.out.html
		json = TALOS.out.json
		panelapp = TALOS.out.panelapp
		labelled_sv = TALOS.out.labelled_sv
		sv_annotated = ch_sv_fresh
}

output {
	// proposed input TSV for the next cycle - published to the root of the output directory
	next_input_tsv {
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
	labelled_sv {
		path { id, _labelled_sv, _labelled_sv_idx -> "${id}_outputs" }
	}
	sv_annotated {
		path { id, _vcf, _vcf_idx -> "${id}_outputs" }
	}
}
