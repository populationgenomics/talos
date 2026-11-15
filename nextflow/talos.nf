#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { AnnotateMitoVcf } from './modules/talos/AnnotateMitoVcf/main'
include { AnnotateVcfWithFreshClinvar } from './modules/talos/AnnotateVcfWithFreshClinvar/main'
include { ConcatLabelledVcfs } from './modules/talos/ConcatLabelledVcfs/main'
include { UnifiedPanelAppParser } from './modules/talos/UnifiedPanelAppParser/main'
include { RunSvFiltering } from './modules/talos/RunSvFiltering/main'
include { RunSmallFiltering } from './modules/talos/RunSmallFiltering/main'
include { ValidateMOI } from './modules/talos/ValidateMOI/main'
include { HPOFlagging } from './modules/talos/HPOFlagging/main'
include { CreateTalosHTML } from './modules/talos/CreateTalosHTML/main'
include { StartupChecks } from './modules/talos/StartupChecks/main'

workflow TALOS {
	take:
		ch_mane
		ch_gff
		ch_ref_genome
		// one entry per cohort: [cohort, meta], meta being the cohort's row Map as built by TALOS_PIPELINE - pedigree
		// and config are Paths, the optional history, ext_ids, seqr_map, mito and str are a Path or [], and
		// analysis_dir/report_dir are where its results and HTML report are published. has_sv says whether an
		// entry for this cohort will arrive on ch_sv_annotated
		ch_meta
		// one entry per annotated shard: [cohort, vcf]
		ch_shards
		// one entry per cohort with SV data: [cohort, vcf, tbi]
		ch_sv_annotated

    main:
    // existence of these files is necessary for starting the workflow
    // we open them as a channel, and pass the channel through to the method
    ch_hpo_file = channel.fromPath(params.hpo, checkIfExists: true).first()
    ch_gen2phen = channel.fromPath(params.gen2phen, checkIfExists: true).first()
    ch_phenio = channel.fromPath(params.phenio_db, checkIfExists: true).first()

    // current year-month as a String, used to prompt for up to date resource updates
    def current_month = new java.util.Date().format('yyyy-MM')
    def timestamp = new java.util.Date().format('yyyy-MM-dd')

    // check if clinvar and panelapp data exist using the timestamp
    // the "all" VCF carries every decision (Benign included), and is the echtvar annotation source
    // for both the small-variant and mito paths - the zip below is the encoded form of it
    String current_clinvarbitration_pm5 = "${params.processed_annotations}/clinvarbitration_${current_month}.pm5.json"
    String current_clinvarbitration_zip = "${params.processed_annotations}/clinvarbitration_${current_month}.zip"

    if (!file(current_clinvarbitration_pm5).exists() || !file(current_clinvarbitration_zip).exists()) {
        println "ClinvArbitration data for this month (=${current_clinvarbitration_pm5}, ${current_clinvarbitration_zip}) doesn't exist, run the Talos Prep workflow"
        exit 1
    }

    // read in each Clinvar input source as channel
    ch_clinvar_pm5 = channel.fromPath(current_clinvarbitration_pm5, checkIfExists: true).first()
    // pre-encoded by the preparation workflow, once per month - fans out across every shard of every cohort
    ch_clinvar_zip = channel.fromPath(current_clinvarbitration_zip, checkIfExists: true).first()

    String panelapp_path = "${params.processed_annotations}/panelapp_${current_month}.json"

    if (!file(panelapp_path).exists()) {
        println "PanelApp data for this month (${panelapp_path}) doesn't exist, run the Talos Prep workflow"
        exit 1
    }
    ch_panelapp = channel.fromPath(panelapp_path, checkIfExists: true).first()

    // run pre-Talos startup checks, tasting the first shard to arrive for each cohort. Grouping the
    // cohort first would delay the checks until every shard was annotated, which is the opposite of
    // what they are for - the cost is that which shard gets checked can vary between runs, so this
    // (cheap) process is not reliably -resume cacheable
    ch_first_shards = ch_shards
        .unique { it -> it[0] }
        .join(ch_meta)
        .map { cohort, vcf, meta -> tuple(cohort, vcf, meta.pedigree, meta.config) }

    StartupChecks(
        ch_first_shards,
    )

    UnifiedPanelAppParser(
        StartupChecks.out
            .join(ch_meta)
            .map { cohort, check_file, meta -> tuple(cohort, check_file, meta.config, meta.pedigree, meta.analysis_dir) },
    	ch_panelapp,
    	ch_hpo_file,
    )

    // [cohort, panelapp_data, meta] - the per-cohort context every downstream stage draws from. Waiting on
    // UnifiedPanelAppParser also means waiting on StartupChecks, which feeds it
    ch_ctx = UnifiedPanelAppParser.out.join(ch_meta)

    // apply this run's ClinVar to every shard, then filter and label each one independently
    AnnotateVcfWithFreshClinvar(
        ch_shards,
        ch_clinvar_zip,
    )

    // combine (not join) - the per-cohort context is repeated across every shard of that cohort
    RunSmallFiltering(
        AnnotateVcfWithFreshClinvar.out
            .combine(ch_ctx, by: 0)
            .map { cohort, vcf, panelapp_data, meta -> tuple(cohort, vcf, panelapp_data, meta.pedigree, meta.config) },
        ch_mane,
        ch_clinvar_pm5,
    )

    // gather the labelled shards back to one VCF per cohort
    ConcatLabelledVcfs(
        RunSmallFiltering.out
            .groupTuple(by: 0)
            .join(ch_meta)
            .map { cohort, vcfs, indexes, meta -> tuple(cohort, vcfs, indexes, meta.analysis_dir) },
    )

    // filter & label any annotated SV VCFs. ch_sv_annotated only carries cohorts that had SV data, so this
    // inner join naturally restricts the process to those cohorts
    RunSvFiltering(
        ch_sv_annotated
            .join(ch_ctx)
            .map { cohort, sv_vcf, sv_idx, panelapp_data, meta ->
                tuple(cohort, sv_vcf, sv_idx, panelapp_data, meta.pedigree, meta.config, meta.analysis_dir)
            },
        ch_mane,
    )

    // one entry per cohort: the filtered SV VCF, or [] for a cohort with no SV data. Which cohorts have SVs is
    // known up front (meta.has_sv), so the rest are given [] directly and never wait on RunSvFiltering.
    // Deliberately not `join(..., remainder: true)`: when the right-hand channel never emits, Nextflow cannot
    // infer its arity and emits the bare cohort key instead of [cohort, null], silently emptying every downstream join
    ch_sv_resolved = RunSvFiltering.out
        .map { cohort, sv_vcf, _sv_idx -> tuple(cohort, sv_vcf) }
        .mix(ch_meta.filter { row -> !row[1].has_sv }.map { row -> tuple(row[0], []) })

    // surprise! It's Mito data! Same shape as the SV stream: annotated for cohorts that have it, [] for the rest
    AnnotateMitoVcf(
        ch_ctx
            .filter { row -> row[2].mito }
            .map { cohort, panelapp_data, meta ->
                tuple(cohort, meta.mito, panelapp_data, meta.pedigree, meta.config,
                      file(params.mitimpact_zip, checkIfExists: true),
                      file(params.mitotip_zip, checkIfExists: true),
                      file(params.napogee_zip, checkIfExists: true))
            },
        ch_ref_genome,
        ch_gff,
        ch_clinvar_zip,
    )

    ch_mito_resolved = AnnotateMitoVcf.out
        .mix(ch_meta.filter { row -> !row[1].mito }.map { row -> tuple(row[0], []) })

    // Validate MOI of all variants
    ValidateMOI(
        ConcatLabelledVcfs.out
            .join(ch_ctx)
            .join(ch_mito_resolved)
            .join(ch_sv_resolved)
            .map { cohort, labelled_vcf, labelled_vcf_index, panelapp_data, meta, anno_mito, anno_sv ->
                tuple(cohort, labelled_vcf, labelled_vcf_index, anno_sv, anno_mito, meta.str, panelapp_data,
                      meta.pedigree, meta.config, meta.history)
            },
        timestamp,
    )

    // Flag any relevant HPO terms
    HPOFlagging(
        ValidateMOI.out
            .join(ch_ctx)
            .map { cohort, talos_result_json, panelapp_data, meta ->
                tuple(cohort, talos_result_json, panelapp_data, meta.config, meta.analysis_dir)
            },
        ch_gen2phen,
        ch_phenio,
        timestamp,
    )

    // Generate HTML report
    CreateTalosHTML(
        HPOFlagging.out
            .join(ch_ctx)
            .map { cohort, result_json, panelapp_data, meta ->
                tuple(cohort, result_json, panelapp_data, meta.config, meta.ext_ids, meta.seqr_map, meta.report_dir)
            },
        timestamp,
    )

    emit:
    	json = HPOFlagging.out
    	// [cohort, html]
    	html = CreateTalosHTML.out
}
