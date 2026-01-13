#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { AnnotateClinvarWithBcftools } from './modules/talos/AnnotateClinvarWithBcftools/main'
include { ConvertSpliceVarDb } from './modules/talos/ConvertSpliceVarDb/main'
include { CreateTalosHTML } from './modules/talos/CreateTalosHTML/main'
include { DownloadClinVarFiles } from './modules/talos/DownloadClinVarFiles/main'
include { DownloadPanelApp } from './modules/talos/DownloadPanelApp/main'
include { HPOFlagging } from './modules/talos/HPOFlagging/main'
include { MakeClinvarbitrationPm5 } from './modules/talos/MakeClinvarbitrationPm5/main'
include { ResummariseRawSubmissions } from './modules/talos/ResummariseRawSubmissions/main'
include { RunHailFiltering } from './modules/talos/RunHailFiltering/main'
include { StartupChecks } from './modules/talos/StartupChecks/main'
include { UnifiedPanelAppParser } from './modules/talos/UnifiedPanelAppParser/main'
include { ValidateMOI } from './modules/talos/ValidateMOI/main'

workflow {
    // existence of these files is necessary for starting the workflow
    // we open them as a channel, and pass the channel through to the method
    // pedigree_channel = Channel.fromPath(params.pedigree)
    ch_hpo_file = Channel.fromPath(params.hpo, checkIfExists: true)
    ch_runtime_config = Channel.fromPath(params.runtime_config, checkIfExists: true)
    ch_gen2phen = Channel.fromPath(params.gen2phen, checkIfExists: true)
    ch_phenio = Channel.fromPath(params.phenio_db, checkIfExists: true)
    ch_mane = Channel.fromPath(params.parsed_mane, checkIfExists: true)
    ch_pedigree = Channel.fromPath(params.pedigree, checkIfExists: true)
    ch_mt = Channel.fromPath(params.matrix_table, checkIfExists: true)
    ch_opt_ids = Channel.fromPath(params.ext_id_map, checkIfExists: true)
    ch_seqr_ids = Channel.fromPath(params.seqr_lookup, checkIfExists: true)

    // may not exist on the first run, will be populated using a dummy file
    ch_previous_results = Channel.fromPath(params.previous_results, checkIfExists: true)

    // does this month's clinvarbitration data exist?
    def current_month = new java.util.Date().format('yyyy-MM')
    String current_clinvarbitration_all = "${params.processed_annotations}/clinvarbitration_${current_month}.ht"
    String current_clinvarbitration_pm5 = "${params.processed_annotations}/clinvarbitration_${current_month}.pm5.ht"

    if (file(current_clinvarbitration_pm5).exists()) {
        ch_clinvar_all = Channel.fromPath(current_clinvarbitration_all)
        ch_clinvar_pm5 = Channel.fromPath(current_clinvarbitration_pm5)
    } else {
        // new workflow elements to go and create it
        DownloadClinVarFiles()

        // reinterpret the results using altered heuristics
        ResummariseRawSubmissions(
            DownloadClinVarFiles.out.variants,
            DownloadClinVarFiles.out.submissions,
        )

        ch_gff = Channel.fromPath(params.ensembl_gff, checkIfExists: true)
        ch_ref_fa = Channel.fromPath(
            params.ref_genome,
            checkIfExists: true,
        )

        // annotate the SNV VCF using BCFtools
        AnnotateClinvarWithBcftools(
            ResummariseRawSubmissions.out.vcf,
            ch_ref_fa,
            ch_gff,
        )

        MakeClinvarbitrationPm5(
            AnnotateClinvarWithBcftools.out.tsv,
        )

        ch_clinvar_all = ResummariseRawSubmissions.out.ht
        ch_clinvar_pm5 = MakeClinvarbitrationPm5.out.ht
    }

    // run pre-Talos startup checks
    StartupChecks(
        ch_mt,
        ch_pedigree,
        ch_clinvar_all,
        ch_runtime_config,
    )

    // download everything in PanelApp - unless it exists from a previous download
    String panelapp_dump = "${params.processed_annotations}/panelapp_${current_month}.json"
    if(file(panelapp_dump).exists()) {
		ch_panelapp = Channel.fromPath(panelapp_dump)
	} else {
		DownloadPanelApp(
			ch_mane,
			ch_runtime_config,
		)
		ch_panelapp = DownloadPanelApp.out.json
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
