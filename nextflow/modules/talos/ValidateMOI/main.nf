process ValidateMOI {
    container params.container

    input:
        tuple val(cohort), path(labelled_vcf), path(labelled_vcf_index), path(sv), path(mito), path(panelapp), path(pedigree), path(talos_config), path(previous_results)
        val timestamp

    output:
        tuple val(cohort), path("${cohort}_results_${timestamp}.json")

	script:
		def history_arg = previous_results.name != 'NO_HISTORY' ? "--previous $previous_results" : ''
		def mito_arg = mito.name != 'NO_MITO' ? "--labelled_mito $mito" : ''
		def mito_idx = mito.name != 'NO_MITO' ? "tabix $mito" : ''
		def sv_arg = sv.name != 'NO_SV' ? "--labelled_sv $sv" : ''
		def sv_idx = sv.name != 'NO_SV' ? "tabix $sv" : ''

        """
        set -euo pipefail

        export TALOS_CONFIG=${talos_config}

        ${mito_idx}
        ${sv_idx}

        python -m talos.validate_moi \
            --labelled_vcf ${labelled_vcf} \
            --panelapp ${panelapp} \
            --pedigree ${pedigree} \
            --output ${cohort}_results_${timestamp}.json $history_arg $mito_arg $sv_arg
        """
}
