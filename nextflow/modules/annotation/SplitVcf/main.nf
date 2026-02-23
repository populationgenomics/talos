process SplitVcf {
    container params.container

    input:
        tuple val(cohort), path(vcf)

    publishDir "${params.outdir}/${cohort}_outputs"

    output:
        tuple val(cohort), path("split_*.bgz")

    script:
    """
    # save the header
    bcftools view -h ${vcf} > header.txt

    # extract the data, split every `vcf_split_n` lines, pipe into a loop that prepends the header and re-compresses.
    bcftools view -H ${vcf} | split -l ${params.vcf_split_n} -d - --filter='cat header.txt - | bgzip -c > split_\${FILE}.vcf.bgz'
    """
}
