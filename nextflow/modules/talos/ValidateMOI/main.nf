process ValidateMOI {
    container params.container

    input:
        tuple val(cohort), path(labelled_vcf), path(labelled_vcf_index), path(sv), path(mito), path(str_vcf), path(panelapp), path(pedigree), path(talos_config), path(previous_results)
        val timestamp

    output:
        tuple val(cohort), path("${cohort}_results_${timestamp}.json")

	// sv, mito, str_vcf and previous_results are optional - absent inputs arrive as [] (falsy, stages nothing)
	script:
		def history_arg = previous_results ? "--previous $previous_results" : ''
		def mito_arg = mito ? "--labelled_mito $mito" : ''
		def mito_idx = mito ? "tabix $mito" : ''
		def sv_arg = sv ? "--labelled_sv $sv" : ''
		def sv_idx = sv ? "tabix $sv" : ''
		def str_arg = str_vcf ? "--str $str_vcf" : ''
		def str_idx = str_vcf ? "tabix -f $str_vcf" : ''

        """
        set -euo pipefail

        export TALOS_CONFIG=${talos_config}

        ${mito_idx}
        ${sv_idx}
        ${str_idx}

        python -m talos.validate_moi \
            --labelled_vcf ${labelled_vcf} \
            --panelapp ${panelapp} \
            --pedigree ${pedigree} \
            --output ${cohort}_results_${timestamp}.json $history_arg $mito_arg $sv_arg $str_arg
        """
}
