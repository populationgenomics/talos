
/*
 * Pipeline parameters - this sets the defaults
 */

params.checkpoint = "$projectDir/assets/NO_FILE"

params.greeting = ["Bonjour", "le", "monde!"]

/*
 * Pipeline parameters
 */
params.input_file = "data/greetings.txt"


process GeneratePanelData {
    // takes the HPO-embellished pedigree and matches panels to participants
    publishDir 'results', mode: 'copy'
    input:
        // the pedigree file
        path pedigree
        // the HPO obo/ontology file
        path hpo

    output:
        path "${params.cohort}_hpo_panel_data.json"

    """
    GeneratePanelData -i ${pedigree} --hpo ${hpo} --out_path ${params.cohort}_hpo_panel_data.json
    """
}

process QueryPanelapp {
    // uses matched panels to query the PanelApp API
    publishDir 'results', mode: 'copy'

    input:
        path hpo_panel_matches

    output:
        path "${params.cohort}_panelapp_results.json"

    // the command
    """
    QueryPanelapp --panels ${panel_data}  --out_path ${params.cohort}_panelapp_results.json
    """
}

process FindGeneSymbolMap {
    // finds corresponding gene symbols for the panelapp results
    publishDir 'results', mode: 'copy'

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
    RunHailFilteringSV \
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
---
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

    GeneratePanelData(params.pedigree, params.hpo)
    // create a channel for inputs from a file
    greeting_ch = Channel.fromPath(params.input_file).splitText() { it.trim() }

    sayHello(greeting_ch)

    // convert the greeting to uppercase
    convertToUpper(sayHello.out)
}