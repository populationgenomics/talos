
process RunHailFiltering {
    container params.container


    // runs the hail small-variant filtering
    publishDir 'results', mode: 'copy'

    input:
        path matrix_table
        path panelapp_data
        path pedigree
        path clinvar
        val checkpoint
        path talos_config

    output:
        tuple \
            path("${params.cohort}_small_variants_labelled.vcf.bgz"), \
            path("${params.cohort}_small_variants_labelled.vcf.bgz.tbi")

    // only write a checkpoint if we were given a path
    script:
        def checkpoint = checkpoint != 'NO_FILE' ? "--checkpoint ${checkpoint}" : ''

    // unzip the ClinvArbitration data directory and MatrixTable
    """
    export TALOS_CONFIG=${talos_config}

    tar -zxf ${clinvar}
    tar -zxf ${matrix_table}

    RunHailFiltering \
        --input ${params.cohort}_small_variants.mt \
        --panelapp ${panelapp_data} \
        --pedigree ${pedigree} \
        --output ${params.cohort}_small_variants_labelled.vcf.bgz \
        --clinvar clinvarbitration_data/clinvar_decisions.ht \
        --pm5 clinvarbitration_data/clinvar_pm5.ht \
        ${checkpoint}
    """
}