#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// deactivated for now
include { ConvertSpliceVarDb } from './modules/talos/ConvertSpliceVarDb/main'

include { ConvertPedToPhenopackets } from './modules/talos/ConvertPedToPhenopackets/main'
include { MakePhenopackets } from './modules/talos/MakePhenopackets/main'
include { GeneratePanelData } from './modules/talos/GeneratePanelData/main'
include { QueryPanelapp } from './modules/talos/QueryPanelapp/main'
include { RunHailFiltering } from './modules/talos/RunHailFiltering/main'
include { ValidateMOI } from './modules/talos/ValidateMOI/main'
include { HPOFlagging } from './modules/talos/HPOFlagging/main'
include { CreateTalosHTML } from './modules/talos/CreateTalosHTML/main'

workflow {
    // existence of these files is necessary for starting the workflow
    // we open them as a channel, and pass the channel through to the method
    // pedigree_channel = Channel.fromPath(params.pedigree)
    ch_hpo_pedigree = channel.fromPath(params.hpo_pedigree, checkIfExists: true)
    ch_hpo_file = channel.fromPath(params.hpo, checkIfExists: true)
    ch_runtime_config = channel.fromPath(params.runtime_config, checkIfExists: true)
    ch_clinvar_tar = channel.fromPath(params.clinvar, checkIfExists: true)
    ch_gen2phen = channel.fromPath(params.gen2phen, checkIfExists: true)
    ch_phenio_gz = channel.fromPath(params.phenio_db, checkIfExists: true)
    ch_mane = channel.fromPath(params.parsed_mane, checkIfExists: true)

    // make a phenopackets file (CPG-specific)
    // commenting this call out as authenticating inside the container is more effort than it's worth right now
    // the validation cohort doesn't have any phenotypes, so we gain nothing from this
    // MakePhenopackets(params.cohort, params.sequencing_type, ch_hpo_file, params.sequencing_tech)

    // instead... make a phenopackets file from pedigree (CPG-specific)
    ConvertPedToPhenopackets(
        ch_hpo_pedigree
    )

    // we can do this by saving as an object, or inside the method call
    GeneratePanelData(
        ConvertPedToPhenopackets.out[1],
        ch_hpo_file
    )

    QueryPanelapp(
        GeneratePanelData.out,
        ch_runtime_config
    )

    // run the hail filtering, using a Tarball'd MT path provided in config
    ch_mt_tar = channel.fromPath(params.matrix_tar, checkIfExists: true)
    RunHailFiltering(
        ch_mt_tar,
        QueryPanelapp.out,
        ConvertPedToPhenopackets.out[0],
        ch_clinvar_tar,
        ch_runtime_config,
    )

    // Validate MOI of all variants
    ValidateMOI(
        RunHailFiltering.out,
        QueryPanelapp.out,
        ConvertPedToPhenopackets.out[0],
        GeneratePanelData.out,
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
        QueryPanelapp.out,
        ch_runtime_config,
    )
}
