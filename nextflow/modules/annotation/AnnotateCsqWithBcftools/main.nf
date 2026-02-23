process AnnotateCsqWithBcftools {
    container params.container

    input:
        tuple val(cohort), path(vcf)
        path gff3
        path reference

    // annotate this VCF with gnomAD data
    publishDir "${params.outdir}/${cohort}_outputs", mode: 'copy'

    output:
        tuple val(cohort), path("${vcf.simpleName}_csq.vcf.bgz")

    script:
    """
    bcftools index -t ${vcf}
    bcftools csq --force -f "${reference}" \
        --local-csq \
        -g ${gff3} \
        --unify-chr-names 'chr,-,chr' \
        -B 20 \
        -Oz -o "${vcf.simpleName}_csq.vcf.bgz" \
        ${vcf}
    """
}
