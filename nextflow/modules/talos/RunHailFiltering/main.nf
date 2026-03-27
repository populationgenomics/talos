process RunHailFiltering {
    container params.container

    input:
        tuple val(cohort), path(mts), path(panelapp_data), path(check_file), path(pedigree), path(talos_config)
        path clinvar_all
        path clinvar_pm5

    output:
        tuple val(cohort), path("${cohort}_small_variants_labelled.vcf.bgz"), path("${cohort}_small_variants_labelled.vcf.bgz.tbi")

    script:
        def mt_string = (mts.collect().size() > 1) ? mts.sort{ it.name } : mts
        """
        set -euo pipefail
        export TALOS_CONFIG=${talos_config}

        python -m talos.run_hail_filtering \
            --input ${mt_string} \
            --panelapp ${panelapp_data} \
            --pedigree ${pedigree} \
            --output ${cohort}_small_variants_labelled.vcf.bgz \
            --clinvar ${clinvar_all} \
            --pm5 ${clinvar_pm5} \
            --checkpoint checkpoint
        """
}
