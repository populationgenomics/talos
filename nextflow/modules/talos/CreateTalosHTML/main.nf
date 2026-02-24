process CreateTalosHTML {
    container params.container

    // generate the HTML report
    publishDir "${params.outdir}/${cohort}_outputs", mode: 'copy'

    input:
        tuple val(cohort), path(talos_result_json), path(panelapp_data), path(talos_config), path(ext_ids), path(seqr_ids)

    def timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm')

    output:
        tuple val(cohort), path("${cohort}_report_${timestamp}.html")

    // add the external IDs file if provided
    script:
        def ext_id_arg = ext_ids.name != 'NO_FILE' ? "--ext_ids $ext_ids" : ''
        def seqr_arg = seqr_ids.name != 'NO_SEQR_FILE' ? "--seqr_ids $seqr_ids" : ''

        """
        export TALOS_CONFIG=${talos_config}
        mkdir ${cohort}_report
        python -m talos.create_talos_html \
            --input ${talos_result_json} \
            --panelapp ${panelapp_data} \
            --output ${cohort}_report_${timestamp}.html $ext_id_arg $seqr_arg
        """
}
