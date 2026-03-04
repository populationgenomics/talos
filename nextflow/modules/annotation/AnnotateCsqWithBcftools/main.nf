process AnnotateCsqWithBcftools {
    container params.container

    input:
        tuple val(cohort), path(vcf)
        path gff3
        path reference

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
