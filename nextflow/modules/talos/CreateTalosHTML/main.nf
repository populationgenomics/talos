
process CreateTalosHTML {
    container params.container

    // generate the HTML report
    publishDir params.output_dir, mode: 'copy'

    input:
        path talos_result_json
        path panelapp_data
        path talos_config
        path ext_ids

    // tarchive of the index page and all per-family reports
    output:
        path "${params.cohort}_report.tar.gz"

	// add the external IDs file if provided
	script:
	def ext_id_arg = ext_ids.name != 'NO_FILE' ? "--ext_ids $ext_ids" : ''

    """
    export TALOS_CONFIG=${talos_config}
    mkdir ${params.cohort}_report
    CreateTalosHTML \
        --input ${talos_result_json} \
        --panelapp ${panelapp_data} \
        --output ${params.cohort}_report/${params.cohort}_report.html $ext_id_arg
    tar --no-xattrs -czf ${params.cohort}_report.tar.gz ${params.cohort}_report
    """
}
