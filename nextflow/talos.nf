#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// deactivated for now
include { ConvertSpliceVarDb } from './modules/talos/ConvertSpliceVarDb/main'
include { ConvertPedToPhenopackets } from './modules/talos/ConvertPedToPhenopackets/main'
include { MakePhenopackets } from './modules/talos/MakePhenopackets/main'
include { DownloadPanelApp } from './modules/talos/DownloadPanelApp/main'
include { UnifiedPanelAppParser } from './modules/talos/UnifiedPanelAppParser/main'
include { RunHailFiltering } from './modules/talos/RunHailFiltering/main'
include { ValidateMOI } from './modules/talos/ValidateMOI/main'
include { HPOFlagging } from './modules/talos/HPOFlagging/main'
include { CreateTalosHTML } from './modules/talos/CreateTalosHTML/main'

workflow {
    // existence of these files is necessary for starting the workflow
    // we open them as a channel, and pass the channel through to the method
    // pedigree_channel = channel.fromPath(params.pedigree)
    ch_hpo_file = channel.fromPath(params.hpo, checkIfExists: true)
    ch_runtime_config = channel.fromPath(params.runtime_config, checkIfExists: true)
    ch_clinvar_tar = channel.fromPath(params.clinvar, checkIfExists: true)
    ch_gen2phen = channel.fromPath(params.gen2phen, checkIfExists: true)
    ch_phenio_gz = channel.fromPath(params.phenio_db, checkIfExists: true)
    ch_mane = channel.fromPath(params.parsed_mane, checkIfExists: true)

    // if phenopackets file is available, read it in
    if (file(params.phenopackets).exists() && file(params.pedigree).exists()) {
    	ch_phenopackets = channel.fromPath(params.phenopackets)
		ch_pedigree = channel.fromPath(params.pedigree)
    }

    // otherwise check if there's a pedigree with HPOs - use that
    else if(file(params.hpo_pedigree).exists()) {
		ch_hpo_pedigree = channel.fromPath(params.hpo_pedigree)
		ConvertPedToPhenopackets(
			ch_hpo_pedigree,
		)
		ch_pedigree = ConvertPedToPhenopackets.out.ped
		ch_phenopackets = ConvertPedToPhenopackets.out.phenopackets
    }

    // fallback, read the plain pedigree in, which will generate a plain phenopacket file
    else {
    	ch_maybe_pedigree = channel.fromPath(
			params.pedigree,
			checkIfExists: true,
		)
		ConvertPedToPhenopackets(
			ch_maybe_pedigree,
		)
		ch_pedigree = ConvertPedToPhenopackets.out.ped
		ch_phenopackets = ConvertPedToPhenopackets.out.phenopackets
    }

    // download everything in PanelApp - unless it exists from a previous download
    if(file(params.panelapp).exists()) {
		ch_panelapp = channel.fromPath(params.panelapp)
	}
	else {
		DownloadPanelApp(
			ch_mane
		)
		ch_panelapp = DownloadPanelApp.out
	}

    // UnifiedPanelAppParser
    UnifiedPanelAppParser(
        ch_runtime_config,
    	ch_panelapp,
    	ch_phenopackets,
    	ch_hpo_file
    )

    // run the hail filtering, using a Tarball'd MT path provided in config
    ch_mt_tar = channel.fromPath(params.matrix_tar, checkIfExists: true)
    RunHailFiltering(
        ch_mt_tar,
        UnifiedPanelAppParser.out,
        ch_pedigree,
        ch_clinvar_tar,
        ch_runtime_config,
    )

    // Validate MOI of all variants
    ValidateMOI(
        RunHailFiltering.out,
        UnifiedPanelAppParser.out,
        ch_pedigree,
        ch_runtime_config,
    )

    // Flag any relevant HPO terms
    HPOFlagging(
        ValidateMOI.out,
        ch_mane,
        ch_gen2phen,
        ch_phenio_gz,
        ch_runtime_config,
    )

    // Generate HTML report - only suited to single-report runs
    CreateTalosHTML(
        HPOFlagging.out.pheno_annotated,
        UnifiedPanelAppParser.out,
        ch_runtime_config,
    )
}
