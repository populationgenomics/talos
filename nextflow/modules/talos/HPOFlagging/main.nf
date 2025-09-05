process HPOFlagging {
    container params.container

    // flag results in the result JSON which are good phenotypic matches
    publishDir params.output_dir, mode: 'copy'

    input:
        path talos_result_json
        path gene_symbol_map
        path gene_to_phenotype
        path phenio_db
        path talos_config
        path previous_results

	def timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm')

    output:
        path "${params.cohort}_full_report_${timestamp}.json"

	script:
		def history_arg = previous_results.name != 'NO_HISTORY' ? "--previous $previous_results" : ''

    """
    export TALOS_CONFIG=${talos_config}
    HPOFlagging \
         --input ${talos_result_json} \
         --mane_json ${gene_symbol_map} \
         --gen2phen ${gene_to_phenotype} \
         --phenio ${phenio_db} \
         --output ${params.cohort}_full_report_${timestamp}.json $history_arg
    """
}
