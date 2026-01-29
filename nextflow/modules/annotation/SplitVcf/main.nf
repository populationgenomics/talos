process SplitVcf {
    container params.container

    input:
        path vcf

    publishDir params.cohort_output_dir

    output:
        path("split_*.bgz")

    script:
    """
    # save the header
    bcftools view -h ${vcf} > header.txt

    # extract the data, split every `vcf_split_n` lines, pipe into a loop that prepends the header and re-compresses.
    bcftools view -H ${vcf} | split -l ${params.vcf_split_n} -d - --filter='cat header.txt - | bgzip -c > split_\${FILE}.vcf.bgz'
    """
}
