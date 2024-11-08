
process FindGeneSymbolMap {
    container params.container

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
