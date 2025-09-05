
process ValidateMOI {
    container params.container


    // process the labelled variants
    publishDir params.output_dir, mode: 'copy'

    input:
        tuple path(labelled_vcf), path(labelled_vcf_index)
        path panelapp
        path pedigree
        path talos_config
        path previous_results

	def timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm')

    output:
        path"${params.cohort}_results_${timestamp}.json"

	script:
		def history_arg = previous_results.name != 'NO_HISTORY' ? "--previous $previous_results" : ''
    """
    export TALOS_CONFIG=${talos_config}

    ValidateMOI \
        --labelled_vcf ${labelled_vcf} \
        --panelapp ${panelapp} \
        --pedigree ${pedigree} \
        --output ${params.cohort}_results_${timestamp}.json $history_arg
    """
}
