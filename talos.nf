
/*
 * Pipeline parameters - this sets the defaults
 */

params.checkpoint = "$projectDir/assets/NO_FILE"

params.greeting = ["Bonjour", "le", "monde!"]

/*
 * Pipeline parameters
 */
params.input_file = "data/greetings.txt"


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

    // only write a checkpoint if we were given a path
    script:
        def checkpoint = checkpoint.name != 'NO_FILE' ? "${checkpoint}" : ''

    """

    RunHailFiltering \
        --mt ${matrix_table} \
        --panelapp ${panelapp_data} \
        --pedigree ${pedigree} \
        --vcf_out ${params.cohort}_small_variants.vcf.bgz \
        --clinvar ${clinvar} \
        --clinvar ${pm5} \
        --checkpoint ${checkpoint}
    """

    output:
        path tuple path("${params.cohort}_small_variants.vcf.bgz"), path("${params.cohort}_small_variants.vcf.bgz.tbi")

    """
    RunHailFilteringSV \
        --mt ${matrix_table} \
        --panelapp ${panelapp_data} \
        --pedigree ${pedigree} \
        --vcf_out ${params.cohort}_small_variants.vcf.bgz
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
        path tuple path("${params.cohort}_small_variants.vcf.bgz"), path("${params.cohort}_small_variants.vcf.bgz.tbi")

    """
    RunHailFilteringSV \
        --mt ${matrix_table} \
        --panelapp ${panelapp_data} \
        --pedigree ${pedigree} \
        --vcf_out ${params.cohort}_small_variants.vcf.bgz
    """
}

process sayHello {
    // this takes the resulting file from the work/hash/hash directory, writes into "results"
    publishDir 'results', mode: 'copy'

    // if doubled
    output:
        tuple path(input_bam), path("${input_bam}.bai")

    // if passed
    input:
        tuple path(input_bam), path(input_bam_index)

    output:
        path "${greeting}-output.txt"

    input:
        val greeting

    """
    echo '$greeting' > '$greeting-output.txt'
    """
}

/*
 * Use a text replace utility to convert the greeting to uppercase
 */
process convertToUpper {

    publishDir 'results', mode: 'copy'

    input:
        path input_file

    output:
        path "UPPER-${input_file}"

    """
    cat '$input_file' | tr '[a-z]' '[A-Z]' > UPPER-${input_file}
    """
}

workflow {
    // existence of these files is necessary for starting the workflow
    // we open them as a channel, and pass the channel through to the method
    pedigree_channel = Channel.fromPath(params.pedigree)
    hpo_file_channel = Channel.fromPath(params.hpo)

    // make a phenopackets file (CPG-specific)
    MakePhenopackets(params.cohort, params.sequencing_type, hpo_file_channel, params.sequencing_tech)

    // we can do this by saving as an object, or inside the method call
    GeneratePanelData(MakePhenopackets.out, hpo_file_channel)

    QueryPanelapp(GeneratePanelData.out, params.runtime_config)
}