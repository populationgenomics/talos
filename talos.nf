#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process VcfToMt {
    publishDir params.output_dir, mode: 'copy'

    input:

        // the path to the HPO file
        path vcf
        path vcf_idx

    output:
        path "${params.cohort}_small_variants.mt.tar.gz"

    """
    VcfToMt \
        --input ${vcf} \
        --output ${params.cohort}_small_variants.mt
    tar -czf ${params.cohort}_small_variants.mt.tar.gz ${params.cohort}_small_variants.mt
    """
}


process ConvertPedToPhenopackets {
    // takes the pedigree file, and converts it to a phenopackets file and regular pedigree file
    publishDir params.output_dir, mode: 'copy'

    input:
        // the pedigree with embedded HPO terms
        path pedigree

    output:
        path "${params.cohort}_pedigree.ped", emit: "ped"
        path "${params.cohort}_phenopackets.json", emit: "phenopackets"

    """
    ConvertPedToPhenopackets \
        --input ${pedigree} \
        --output ${params.cohort}
    """
}


process MakePhenopackets {
    publishDir params.output_dir, mode: 'copy'

    input:
        // the name of the dataset for use in API queries
        val dataset

        // the sequencing type
        val sequencing_type

        // the path to the HPO file
        path hpo

        // sequencing technology, not currently overridden
        val sequencing_tech

    output:
        path "${params.cohort}_phenopackets.json"

    """
    MakePhenopackets \
        --dataset ${dataset} \
        --output ${params.cohort}_phenopackets.json \
        --type ${sequencing_type} \
        --hpo ${hpo} \
        --tech ${sequencing_tech}
    """
}

process GeneratePanelData {
    // takes the HPO-embellished pedigree and matches panels to participants
    publishDir params.output_dir, mode: 'copy'

    input:
        // the pedigree file
        path pedigree
        // the HPO obo/ontology file
        path hpo

    output:
        path "${params.cohort}_hpo_panel_data.json"

    """
    GeneratePanelData --input ${pedigree} --hpo ${hpo} --output ${params.cohort}_hpo_panel_data.json
    """
}

process QueryPanelapp {
    // uses matched panels to query the PanelApp API
    publishDir params.output_dir, mode: 'copy'

    input:
        path hpo_panel_matches
        path talos_config

    output:
        path "${params.cohort}_panelapp_results.json"

    // the command
    """
    export TALOS_CONFIG=${talos_config}
    QueryPanelapp --input ${hpo_panel_matches}  --output ${params.cohort}_panelapp_results.json
    """
}

process FindGeneSymbolMap {
    // finds corresponding gene symbols for the panelapp results
    publishDir params.output_dir, mode: 'copy'

    input:
        path panelapp_data

    output:
        path "${params.cohort}_symbol_to_ensg.json"

    // the command
    """
    FindGeneSymbolMap --panelapp ${panel_data}  --out_path ${params.cohort}_symbol_to_ensg.json
    """
}

process RunHailFiltering {
    // runs the hail small-variant filtering
    publishDir 'results', mode: 'copy'

    input:
        path matrix_table
        path panelapp_data
        path pedigree
        path clinvar
        val checkpoint
        path talos_config

    output:
        tuple \
            path("${params.cohort}_small_variants_labelled.vcf.bgz"), \
            path("${params.cohort}_small_variants_labelled.vcf.bgz.tbi")

    // only write a checkpoint if we were given a path
    script:
        def checkpoint = checkpoint != 'NO_FILE' ? "--checkpoint ${checkpoint}" : ''

    // unzip the ClinvArbitration data directory and MatrixTable
    """
    export TALOS_CONFIG=${talos_config}

    tar -zxf ${clinvar}
    tar -zxf ${matrix_table}

    RunHailFiltering \
        --input ${params.cohort}_small_variants.mt \
        --panelapp ${panelapp_data} \
        --pedigree ${pedigree} \
        --output ${params.cohort}_small_variants_labelled.vcf.bgz \
        --clinvar clinvarbitration_data/clinvar_decisions.ht \
        --pm5 clinvarbitration_data/clinvar_pm5.ht \
        ${checkpoint}
    """
}

process RunHailFilteringSV {
    // runs the hail structural-variant filtering
    // here we can have multiple files - should we derive output name from input?
    publishDir 'results', mode: 'copy'

    input:
        path matrix_table
        path panelapp_data
        path pedigree

    output:
        path "${params.cohort}_small_variants.vcf.bgz", emit: "vcf"
        path "${params.cohort}_small_variants.vcf.bgz.tbi", emit: "index"

    """
    RunHailFilteringSV \
        --mt ${matrix_table} \
        --panelapp ${panelapp_data} \
        --pedigree ${pedigree} \
        --vcf_out ${params.cohort}_small_variants.vcf.bgz
    """
}

process ValidateMOI {
    // process the labelled variants
    publishDir params.output_dir, mode: 'copy'

    input:
        tuple path(labelled_vcf), path(labelled_vcf_index)
        path panelapp
        path pedigree
        path hpo_panel_matches
        path talos_config

    output:
        path"${params.cohort}_results.json"

    """
    export TALOS_CONFIG=${talos_config}

    ValidateMOI \
        --labelled_vcf ${labelled_vcf} \
        --panelapp ${panelapp} \
        --pedigree ${pedigree} \
        --participant_panels ${hpo_panel_matches} \
        --output ${params.cohort}_results.json
    """
}

workflow {
    // existence of these files is necessary for starting the workflow
    // we open them as a channel, and pass the channel through to the method
    // pedigree_channel = Channel.fromPath(params.pedigree)
    hpo_pedigree_channel = Channel.fromPath(params.hpo_pedigree)
    hpo_file_channel = Channel.fromPath(params.hpo)
    runtime_config_channel = Channel.fromPath(params.runtime_config)
    clinvar_tar_channel = Channel.fromPath(params.clinvar)

    // turn the VCF into a MatrixTable
    input_vcf = Channel.fromPath(params.annotated_vcf)
    input_vcf_idx = Channel.fromPath("${params.annotated_vcf}.tbi")
    VcfToMt(input_vcf, input_vcf_idx)

    // make a phenopackets file from pedigree (CPG-specific)
    ConvertPedToPhenopackets(hpo_pedigree_channel)

    // make a phenopackets file (CPG-specific)
    MakePhenopackets(params.cohort, params.sequencing_type, hpo_file_channel, params.sequencing_tech)

    // we can do this by saving as an object, or inside the method call
    GeneratePanelData(ConvertPedToPhenopackets.out[1], hpo_file_channel)

    QueryPanelapp(
        GeneratePanelData.out,
        runtime_config_channel,
    )

    // run the hail filtering
    RunHailFiltering(
        VcfToMt.out,
        QueryPanelapp.out,
        ConvertPedToPhenopackets.out[0],
        clinvar_tar_channel,
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

}
