#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ConvertSpliceVarDb } from './modules/talos/ConvertSpliceVarDb/main'
include { ConvertPedToPhenopackets } from './modules/talos/ConvertPedToPhenopackets/main'
include { MakePhenopackets } from './modules/talos/MakePhenopackets/main'
include { GeneratePanelData } from './modules/talos/GeneratePanelData/main'
include { QueryPanelapp } from './modules/talos/QueryPanelapp/main'
include { FindGeneSymbolMap } from './modules/talos/FindGeneSymbolMap/main'
include { RunHailFiltering } from './modules/talos/RunHailFiltering/main'
include { ValidateMOI } from './modules/talos/ValidateMOI/main'
include { HPOFlagging } from './modules/talos/HPOFlagging/main'
include { CreateTalosHTML } from './modules/talos/CreateTalosHTML/main'

workflow {
    // existence of these files is necessary for starting the workflow
    // we open them as a channel, and pass the channel through to the method
    // pedigree_channel = Channel.fromPath(params.pedigree)
    hpo_pedigree_channel = Channel.fromPath(params.hpo_pedigree)
    hpo_file_channel = Channel.fromPath(params.hpo)
    runtime_config_channel = Channel.fromPath(params.runtime_config)
    clinvar_tar_channel = Channel.fromPath(params.clinvar)
    exomiser_tar_channel = Channel.fromPath(params.exomiser)
    gen2phen_channel = Channel.fromPath(params.gen2phen)
    phenio_db_channel = Channel.fromPath(params.phenio_db)
    svdb_tsv_channel = Channel.fromPath(params.svdb_tsv)

    // convert the SVDB TSV into a Hail Table
    ConvertSpliceVarDb(svdb_tsv_channel)

    // TODO turn the VCF into a MatrixTable
//     input_vcf = Channel.fromPath(params.annotated_vcf).map{ it -> [file(it), file("${it}.tbi")]}
//     VcfToMt(input_vcf)

    // make a phenopackets file from pedigree (CPG-specific)
    ConvertPedToPhenopackets(hpo_pedigree_channel)

    // make a phenopackets file (CPG-specific)
    // commenting this call out as authenticating inside the container is more effort than it's worth right now
    // the validation cohort doesn't have any phenotypes, so we gain nothing from this
    // MakePhenopackets(params.cohort, params.sequencing_type, hpo_file_channel, params.sequencing_tech)

    // we can do this by saving as an object, or inside the method call
    GeneratePanelData(ConvertPedToPhenopackets.out[1], hpo_file_channel)

    QueryPanelapp(GeneratePanelData.out, runtime_config_channel)

    FindGeneSymbolMap(QueryPanelapp.out, runtime_config_channel)

    // run the hail filtering
    RunHailFiltering(
        VcfToMt.out,
        QueryPanelapp.out,
        ConvertPedToPhenopackets.out[0],
        clinvar_tar_channel,
        exomiser_tar_channel,
        ConvertSpliceVarDb.out,
        params.checkpoint,
        runtime_config_channel,
    )

    // Validate MOI of all variants
    ValidateMOI(
        RunHailFiltering.out,
        QueryPanelapp.out,
        ConvertPedToPhenopackets.out[0],
        GeneratePanelData.out,
        runtime_config_channel,
    )

    // Flag any relevant HPO terms
    HPOFlagging(
        ValidateMOI.out,
        FindGeneSymbolMap.out,
        gen2phen_channel,
        phenio_db_channel,
        runtime_config_channel,
    )

    // Generate HTML report - only suited to single-report runs
    CreateTalosHTML(
        HPOFlagging.out.pheno_annotated,
        QueryPanelapp.out,
        runtime_config_channel,
    )
}
