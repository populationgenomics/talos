
process CreateTalosHTML {
    container params.container

    // generate the HTML report
    publishDir params.output_dir, mode: 'copy'

    input:
        path talos_result_json
        path panelapp_data
        path talos_config

    output:
        path "${params.cohort}_report.html"

    """
    export TALOS_CONFIG=${talos_config}
    CreateTalosHTML \
        --input ${talos_result_json} \
        --panelapp ${panelapp_data} \
        --output ${params.cohort}_report.html
    """
}
