process RunHailFiltering {
    container params.container

    publishDir params.output_dir, mode: 'copy'

    input:
        path mts
        path panelapp_data
        path pedigree
        path clinvar_all
        path clinvar_pm5
        path talos_config
        path check_file

    output:
        tuple \
            path("${params.cohort}_small_variants_labelled.vcf.bgz"), \
            path("${params.cohort}_small_variants_labelled.vcf.bgz.tbi")

    script:
        def mt_string = (mts.collect().size() > 1) ? mts.sort{ it.name } : mts
        """
        export TALOS_CONFIG=${talos_config}

        python -m talos.run_hail_filtering \
            --input ${mt_string} \
            --panelapp ${panelapp_data} \
            --pedigree ${pedigree} \
            --output ${params.cohort}_small_variants_labelled.vcf.bgz \
            --clinvar ${clinvar_all} \
            --pm5 ${clinvar_pm5} \
            --checkpoint checkpoint
        """
}
