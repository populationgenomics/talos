
process VcfToMt {
    publishDir params.output_dir, mode: 'copy'

    input:
        // input VCF, matched with its index
        path vcf

    output:
        path "${params.cohort}_small_variants.mt.tar.gz"

    """
    VcfToMt \
        --input ${vcf[0]} \
        --output ${params.cohort}_small_variants.mt
    tar -czf ${params.cohort}_small_variants.mt.tar.gz ${params.cohort}_small_variants.mt
    """
}