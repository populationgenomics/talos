
process RunHailFiltering {
    container params.container

    // runs the hail small-variant filtering
    publishDir params.output_dir, mode: 'copy'

    input:
        path mt
        path panelapp_data
        path pedigree
        path clinvar
        path talos_config

    output:
        tuple \
            path("${params.cohort}_small_variants_labelled.vcf.bgz"), \
            path("${params.cohort}_small_variants_labelled.vcf.bgz.tbi")

    // unzip the ClinvArbitration data directory and MatrixTable
    """
    export TALOS_CONFIG=${talos_config}

    tar -zxf ${clinvar}

    RunHailFiltering \
        --input ${mt} \
        --panelapp ${panelapp_data} \
        --pedigree ${pedigree} \
        --output ${params.cohort}_small_variants_labelled.vcf.bgz \
        --clinvar clinvarbitration_data/clinvar_decisions.ht \
        --pm5 clinvarbitration_data/clinvar_pm5.ht \
        --checkpoint checkpoint
    """
}
