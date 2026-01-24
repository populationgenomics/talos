process AnnotateWithEchtvar {
    container params.container

    input:
        tuple path(vcf), path(tbi)
        path gnomad_zip
        path am_zip

    publishDir params.cohort_output_dir

    output:
        path("${vcf.simpleName}_echtvar.vcf.bgz")

    script:
        """
        set -ex
        echtvar anno \
            -e ${gnomad_zip} \
            -e ${am_zip} \
            -i "gnomad_AF_joint < 0.05" \
            ${vcf} \
            "${vcf.simpleName}_echtvar.vcf.bgz"
        """
}
