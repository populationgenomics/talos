
process FilterVcfToBedWithBcftools {
    container params.container

    // take the merged VCF and index
    input:
        path(vcf)
        path(tbi)
        path (bed_file)

    // strip this VCF down to sites-only
    publishDir params.cohort_output_dir, mode: 'copy'

    output:
        tuple \
            path("${params.cohort}_merged_filtered.vcf.bgz"), \
            path("${params.cohort}_merged_filtered.vcf.bgz.tbi")

    script:
    """
    bcftools view \
        --write-index=tbi \
         -R ${bed_file} \
        -Oz -o "${params.cohort}_merged_filtered.vcf.bgz" \
        ${vcf}
    """
}
