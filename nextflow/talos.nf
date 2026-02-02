#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// deactivated for now
include { ConvertSpliceVarDb } from './modules/talos/ConvertSpliceVarDb/main'
include { UnifiedPanelAppParser } from './modules/talos/UnifiedPanelAppParser/main'
include { RunHailFiltering } from './modules/talos/RunHailFiltering/main'
include { ValidateMOI } from './modules/talos/ValidateMOI/main'
include { HPOFlagging } from './modules/talos/HPOFlagging/main'
include { CreateTalosHTML } from './modules/talos/CreateTalosHTML/main'
include { StartupChecks } from './modules/talos/StartupChecks/main'

workflow {
    main :
    // existence of these files is necessary for starting the workflow
    // we open them as a channel, and pass the channel through to the method
    // pedigree_channel = Channel.fromPath(params.pedigree)
    ch_hpo_file = Channel.fromPath(params.hpo, checkIfExists: true)
    ch_runtime_config = Channel.fromPath(params.runtime_config, checkIfExists: true)
    ch_gen2phen = Channel.fromPath(params.gen2phen, checkIfExists: true)
    ch_phenio = Channel.fromPath(params.phenio_db, checkIfExists: true)
    ch_mane = Channel.fromPath(params.parsed_mane, checkIfExists: true)
    ch_pedigree = Channel.fromPath(params.pedigree, checkIfExists: true)
    ch_opt_ids = Channel.fromPath(params.ext_id_map, checkIfExists: true)
    ch_seqr_ids = Channel.fromPath(params.seqr_lookup, checkIfExists: true)

    // find all matrix tables in the cohort output directory - require at least one
    ch_mts = Channel.fromPath("${params.output_dir}/*.mt", type: 'dir', checkIfExists: true).collect()

    // may not exist on the first run, will be populated using a dummy file
    ch_previous_results = Channel.fromPath(params.previous_results, checkIfExists: true)

    // current year-month as a String, used to prompt for up to date resource updates
    def current_month = new java.util.Date().format('yyyy-MM')

    // check if clinvar and panelapp data exist using the timestamp
    String current_clinvarbitration_all = "${params.processed_annotations}/clinvarbitration_${current_month}.ht"
    String current_clinvarbitration_pm5 = "${params.processed_annotations}/clinvarbitration_${current_month}.pm5.ht"

    if (!file(current_clinvarbitration_pm5).exists()) {
        println "ClinvArbitration data for this month (${current_clinvarbitration_pm5}) doesn't exist, run the Talos Prep workflow"
        exit 1
    }

    // read in each Clinvar input source as channel
    ch_clinvar_all = Channel.fromPath(current_clinvarbitration_all, checkIfExists: true)
    ch_clinvar_pm5 = Channel.fromPath(current_clinvarbitration_pm5, checkIfExists: true)

    String panelapp = "${params.processed_annotations}/panelapp_${current_month}.json"

    if (!file(panelapp).exists()) {
        println "PanelApp data for this month (${panelapp}) doesn't exist, run the Talos Prep workflow"
        exit 1
    }
    ch_panelapp = Channel.fromPath(panelapp, checkIfExists: true)

    // run pre-Talos startup checks
    StartupChecks(
        ch_mts,
        ch_pedigree,
        ch_clinvar_all,
        ch_runtime_config,
    )

    // UnifiedPanelAppParser
    UnifiedPanelAppParser(
        ch_runtime_config,
    	ch_panelapp,
    	ch_pedigree,
    	ch_hpo_file,
    	StartupChecks.out,
    )

    RunHailFiltering(
        ch_mts,
        UnifiedPanelAppParser.out,
        ch_pedigree,
        ch_clinvar_all,
        ch_clinvar_pm5,
        ch_runtime_config,
        StartupChecks.out,
    )

    // Validate MOI of all variants
    ValidateMOI(
        RunHailFiltering.out,
        UnifiedPanelAppParser.out,
        ch_pedigree,
        ch_runtime_config,
        ch_previous_results,
    )

    // Flag any relevant HPO terms
    HPOFlagging(
        ValidateMOI.out,
        ch_mane,
        ch_gen2phen,
        ch_phenio,
        ch_runtime_config,
    )

    // Generate HTML report - only suited to single-report runs
    CreateTalosHTML(
        HPOFlagging.out,
        UnifiedPanelAppParser.out,
        ch_runtime_config,
        ch_opt_ids,
        ch_seqr_ids,
    )
}
