#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
    Talos Workflow - CPG entrypoint
    ===============================

    Runs the shared TALOS_PIPELINE (nextflow/talos_pipeline.nf) over a set of Metamist cohorts, each publishing
    to its own dataset bucket:

        gs://cpg-<dataset>-{main|test}/<sequencing_type>/talos/<cohort>/annotated/
        gs://cpg-<dataset>-{main|test}/<sequencing_type>/talos/<cohort>/analysis_YYYYMMDD/
        gs://cpg-<dataset>-{main|test}-web/<sequencing_type>/talos/<cohort>/analysis_YYYYMMDD/   (HTML report only)

    The first process, InputTsvFromMetamist, finds each cohort's inputs in Metamist and generates the rest
    (pedigree, external IDs, seqr map, per-cohort config) into that cohort's analysis_YYYYMMDD/inputs/. Its TSV
    then feeds TALOS_PIPELINE exactly as main.nf's input TSV would. There is no next-cycle TSV: history is the
    dataset's latest aip-results analysis in Metamist. With --metamist, each cohort's results JSON is registered
    there as aip-results (the next run's history) and its HTML report as aip-report - see docs/CpgDeployment.md.

    Outputs are published with publishDir, not workflow outputs - see docs/adr/0001.

    Usage:
    nextflow run cpg_main.nf --cohort_ids COH0001,COH0002 --access_level test --sequencing_type genome \
        --talos_config gs://.../talos.toml [--metamist] [other params...]
*/

include { InputTsvFromMetamist } from './nextflow/modules/cpg/InputTsvFromMetamist/main'
include { RegisterAnalysis as RegisterResults } from './nextflow/modules/cpg/RegisterAnalysis/main'
include { RegisterAnalysis as RegisterReport } from './nextflow/modules/cpg/RegisterAnalysis/main'
include { TALOS_PIPELINE } from './nextflow/talos_pipeline'

// fixed once per run - InputTsvFromMetamist writes each cohort's generated inputs under analysis_<runDate>, so the
// published results must land in the same directory
def runDate() {
    workflow.start.format(java.time.format.DateTimeFormatter.ofPattern('yyyyMMdd'))
}

// a list (e.g. from a params file) or a comma-delimited string, as a list of unique IDs
def cohortIds() {
    def raw = params.cohort_ids instanceof List ? params.cohort_ids : params.cohort_ids.toString().split(',')
    raw.collect { cohort_id -> cohort_id.toString().trim() }.findAll { cohort_id -> cohort_id }.unique()
}

// must agree with cohort_root() in talos.cpg_internal.input_tsv_from_metamist. Plain Strings, never a Path.
// suffix selects a sibling bucket, e.g. '-web' for the HTML report
def cohortRoot(dataset, cohort, suffix = '') {
    def bucket = "gs://cpg-${dataset}-${params.access_level == 'test' ? 'test' : 'main'}${suffix}"
    "${bucket}/${params.sequencing_type}/talos/${cohort}".toString()
}

def annotatedDir(dataset, cohort) {
    "${cohortRoot(dataset, cohort)}/annotated".toString()
}

def analysisDir(dataset, cohort) {
    "${cohortRoot(dataset, cohort)}/analysis_${runDate()}".toString()
}

// the HTML report goes to the dataset's -web bucket, mirroring analysisDir
def reportDir(dataset, cohort) {
    "${cohortRoot(dataset, cohort, '-web')}/analysis_${runDate()}".toString()
}

// the Metamist project a dataset's analyses are registered in
def metamistProject(dataset) {
    params.access_level == 'test' ? "${dataset}-test".toString() : dataset
}

// [cohort, project, published URI, meta] for RegisterAnalysis, from a TALOS_PIPELINE [row, published URI] emission
def registration(row, published) {
    tuple(row.cohort, metamistProject(row.cpg_dataset), published, [sequencing_type: params.sequencing_type])
}


workflow {
	main:
	if (!params.cohort_ids || !cohortIds()) {
		println "Required --cohort_ids argument not provided"
		exit 1
	}

	if (!(params.access_level in ['test', 'standard', 'full'])) {
		println "--access_level must be one of test, standard, full (got ${params.access_level})"
		exit 1
	}

	if (!(params.sequencing_type in ['genome', 'exome'])) {
		println "--sequencing_type must be one of genome, exome (got ${params.sequencing_type})"
		exit 1
	}

	if (!params.talos_config) {
		println "Required --talos_config argument not provided"
		exit 1
	}

	InputTsvFromMetamist(
		cohortIds(),
		params.access_level,
		params.sequencing_type,
		file(params.talos_config, checkIfExists: true),
		runDate(),
	)

	// one Map per cohort: the generated TSV row, plus where that cohort's outputs are published
	ch_cohorts = InputTsvFromMetamist.out
		.splitCsv(header: true, sep: '\t')
		.map { row -> row + [
			annotated_dir: annotatedDir(row.cpg_dataset, row.cohort),
			analysis_dir: analysisDir(row.cpg_dataset, row.cohort),
			report_dir: reportDir(row.cpg_dataset, row.cohort),
		] }

	TALOS_PIPELINE(ch_cohorts)

	// register each cohort's published results JSON (aip-results, the next run's history) and HTML report
	// (aip-report) against its Metamist cohort. Each waits for its own file to be published, independently
	if (params.metamist) {
		RegisterResults(
			TALOS_PIPELINE.out.json.map { row, published -> registration(row, published) },
			'aip-results',
			params.metamist_wait_minutes,
		)

		RegisterReport(
			TALOS_PIPELINE.out.html.map { row, published -> registration(row, published) },
			'aip-report',
			params.metamist_wait_minutes,
		)
	}
}
