#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
    Talos Workflow - generic entrypoint
    ===================================

    Runs the shared TALOS_PIPELINE (nextflow/talos_pipeline.nf) over the cohorts in an input TSV:

    1. ANNOTATION: Annotates VCF file(s) with the prepared data.
    2. SV_ANNOTATION: Annotates a joint-called SV VCF, for each cohort with an `sv` VCF.
    3. TALOS: Runs the core Talos analysis/filtering/reporting.

    Annotation (small-variant and SV alike) is idempotent: cohorts whose annotation products are
    already published (keyed on the shard manifest / annotated SV VCF) are reused, not re-annotated.
    A reanalysis cycle is simply re-running this workflow; delete a cohort's `${cohort}_annotated`
    directory to force re-annotation.

    Outputs are published with publishDir, not workflow outputs - see docs/adr/0001. This entrypoint
    decides each cohort's output directories and passes them down in the cohort tuples; nothing below
    it knows the output layout.

    Usage:
    nextflow run main.nf --input_tsv [path] --outdir [path] [other params...]
*/

include { TALOS_PIPELINE } from './nextflow/talos_pipeline'

// analysis outputs are published to a per-run dated directory, ${cohort}_analysis_YYYYMMDD -
// annotation products go to the undated ${cohort}_annotated, which is reused across cycles
def runDate() {
    workflow.start.format(java.time.format.DateTimeFormatter.ofPattern('yyyyMMdd'))
}

// output directories are carried as plain Strings built from params.outdir (itself a String), never
// by interpolating a Path - that drops the gs:// scheme from cloud paths
def outRoot() {
    params.outdir.toString().replaceAll('/+$', '')
}

def annotatedDir(cohort) {
    "${outRoot()}/${cohort}_annotated".toString()
}

def analysisDir(cohort) {
    "${outRoot()}/${cohort}_analysis_${runDate()}".toString()
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
	if (!params.outdir) {
		println "Required --outdir argument not provided"
		exit 1
	}

	if (!params.input_tsv) {
		println "Required --input_tsv argument not provided"
		exit 1
	}

	// one Map per cohort: the TSV row, plus where that cohort's outputs are published
	ch_cohorts = channel.fromPath(params.input_tsv)
		.splitCsv(header: true, sep: '\t')
		.map { row -> row + [annotated_dir: annotatedDir(row.cohort), analysis_dir: analysisDir(row.cohort)] }

	TALOS_PIPELINE(ch_cohorts)

	// a ready-to-use input TSV for the next reanalysis cycle - identical to this run's input, with
	// `history` repointed at the results JSON published for that cohort by this run, written to the
	// root of the output directory. The paths are the published destinations, which don't exist until
	// the run completes, so this is a proposal to be reviewed, not an output consumed by anything in this run
	TALOS_PIPELINE.out.json
		.map { row, published -> tuple(row.cohort, published) }
		.collect(flat: false)
		.map { pairs -> nextInputTsv(file(params.input_tsv), pairs.collectEntries()) }
		.collectFile(name: "talos_input_${runDate()}.tsv", storeDir: file(outRoot()))
}
