#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// deactivated for now
include { ConvertSpliceVarDb } from './modules/talos/ConvertSpliceVarDb/main'
include { DownloadPanelApp } from './modules/talos/DownloadPanelApp/main'
include { UnifiedPanelAppParser } from './modules/talos/UnifiedPanelAppParser/main'
include { RunHailFiltering } from './modules/talos/RunHailFiltering/main'
include { ValidateMOI } from './modules/talos/ValidateMOI/main'
include { HPOFlagging } from './modules/talos/HPOFlagging/main'
include { CreateTalosHTML } from './modules/talos/CreateTalosHTML/main'
include { StartupChecks } from './modules/talos/StartupChecks/main'

workflow {
    // existence of these files is necessary for starting the workflow
    // we open them as a channel, and pass the channel through to the method
    // pedigree_channel = channel.fromPath(params.pedigree)
    ch_hpo_file = channel.fromPath(params.hpo, checkIfExists: true)
    ch_runtime_config = channel.fromPath(params.runtime_config, checkIfExists: true)
    ch_clinvar_tar = channel.fromPath(params.clinvar, checkIfExists: true)
    ch_gen2phen = channel.fromPath(params.gen2phen, checkIfExists: true)
    ch_phenio = channel.fromPath(params.phenio_db, checkIfExists: true)
    ch_mane = channel.fromPath(params.parsed_mane, checkIfExists: true)
    ch_pedigree = channel.fromPath(params.pedigree, checkIfExists: true)
    ch_mt = channel.fromPath(params.matrix_table, checkIfExists: true)
    ch_opt_ids = channel.fromPath(params.ext_id_map, checkIfExists: true)
    ch_seqr_ids = channel.fromPath(params.seqr_lookup, checkIfExists: true)

    // may not exist on the first run, will be populated using a dummy file
    ch_previous_results = channel.fromPath(params.previous_results, checkIfExists: true)

    // run pre-Talos startup checks
    StartupChecks(
        ch_mt,
        ch_pedigree,
        ch_clinvar_tar,
        ch_runtime_config,
    )

    // download everything in PanelApp - unless it exists from a previous download
    if(file(params.panelapp).exists()) {
		ch_panelapp = channel.fromPath(params.panelapp)
	}
	else {
		DownloadPanelApp(
			ch_mane,
			ch_runtime_config,
		)
		ch_panelapp = DownloadPanelApp.out
	}

    // UnifiedPanelAppParser
    UnifiedPanelAppParser(
        ch_runtime_config,
    	ch_panelapp,
    	ch_pedigree,
    	ch_hpo_file,
    	StartupChecks.out,
    )

    RunHailFiltering(
        ch_mt,
        UnifiedPanelAppParser.out,
        ch_pedigree,
        ch_clinvar_tar,
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
