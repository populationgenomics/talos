
process AnnotateGnomadAfWithEchtvar {
    container params.container

    input:
        tuple path(vcf), path(tbi)
        path(zip)

    // annotate this VCF with gnomAD data
    publishDir params.cohort_output_dir, mode: 'copy'

    output:
        path("${params.cohort}_gnomad.vcf.bgz")

    script:
        """
        set -ex
        echtvar anno \
            -e ${zip} \
            ${vcf} \
            "${params.cohort}_gnomad.vcf.bgz"
        """
}
