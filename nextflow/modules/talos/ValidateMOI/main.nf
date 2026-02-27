
process ValidateMOI {
    container params.container

    publishDir "${params.outdir}/${cohort}_outputs", mode: 'copy'

    input:
        tuple val(cohort), path(labelled_vcf), path(labelled_vcf_index), path(panelapp), path(pedigree), path(talos_config), path(previous_results)
        val timestamp

    output:
        tuple val(cohort), path("${cohort}_results_${timestamp}.json")

	script:
		def history_arg = previous_results.name != 'NO_HISTORY' ? "--previous $previous_results" : ''

    """
    export TALOS_CONFIG=${talos_config}

    python -m talos.validate_moi \
        --labelled_vcf ${labelled_vcf} \
        --panelapp ${panelapp} \
        --pedigree ${pedigree} \
        --output ${cohort}_results_${timestamp}.json $history_arg
    """
}
