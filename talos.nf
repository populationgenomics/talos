
/*
 * Pipeline parameters - this sets the defaults
 */

params.checkpoint = "$projectDir/assets/NO_FILE"

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
        env TALOS_CONFIG

    output:
        path "${params.cohort}_panelapp_results.json"

    // the command
    """
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
        path pm5
        path checkpoint
        env TALOS_CONFIG

//     output:
//         path '${params.cohort}_small_variants_labelled.vcf.bgz', emit: 'vcf'
//         path '${params.cohort}_small_variants_labelled.vcf.bgz.tbi', emit: 'index'

    output:
        tuple \
            path("${params.cohort}_small_variants_labelled.vcf.bgz"), \
            path("${params.cohort}_small_variants_labelled.vcf.bgz.tbi")

    // only write a checkpoint if we were given a path
    script:
        def checkpoint = checkpoint.name != 'NO_FILE' ? "--checkpoint ${checkpoint}" : ''

    """
    # unzip the ClinvArbitration data directories
    tar -zxf ${clinvar}
    tar -zxf ${pm5}
    tar -zxf ${matrix_table}

    RunHailFiltering \
        --input ${params.cohort}_small_variants.mt \
        --panelapp ${panelapp_data} \
        --pedigree ${pedigree} \
        --output ${params.cohort}_small_variants_labelled.vcf.bgz \
        --clinvar clinvar_decisions_fake.ht \
        --pm5 clinvar_pm5_fake.ht \
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
        path '${params.cohort}_small_variants.vcf.bgz', emit: 'vcf'
        path '${params.cohort}_small_variants.vcf.bgz.tbi', emit: 'index'

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
        env TALOS_CONFIG

    output:
        path"${params.cohort}_results.json"

    """
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
    pedigree_channel = Channel.fromPath(params.pedigree)
    hpo_file_channel = Channel.fromPath(params.hpo)

    // turn the VCF into a MatrixTable
    input_vcf = Channel.fromPath(params.annotated_vcf)
    input_vcf_idx = Channel.fromPath("${params.annotated_vcf}.tbi")
    VcfToMt(input_vcf, input_vcf_idx)

    // make a phenopackets file (CPG-specific)
    MakePhenopackets(params.cohort, params.sequencing_type, hpo_file_channel, params.sequencing_tech)

    // we can do this by saving as an object, or inside the method call
    GeneratePanelData(MakePhenopackets.out, hpo_file_channel)

    QueryPanelapp(GeneratePanelData.out, params.runtime_config)

    // run the hail filtering
    RunHailFiltering(
        VcfToMt.out,
        QueryPanelapp.out,
        params.pedigree,
        params.clinvar_decisions,
        params.clinvar_pm5,
        params.checkpoint,
        params.runtime_config,
    )

    // Validate MOI of all variants
    ValidateMOI(
        RunHailFiltering.out,
        QueryPanelapp.out,
        params.pedigree,
        GeneratePanelData.out,
        params.runtime_config,
    )

}