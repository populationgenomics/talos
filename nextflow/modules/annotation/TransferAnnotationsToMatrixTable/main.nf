
process TransferAnnotationsToMatrixTable {
    container params.container

    input:
        path(compressed_ht)
        tuple path(vcf), path(tbi)

    publishDir params.cohort_output_dir, mode: 'copy'

    output:
        path("${params.cohort}.mt.tar.gz")

    script:
        """
        set -ex
        tar -xf ${compressed_ht}
        TransferAnnotationsToMatrixTable \
            --input ${vcf} \
            --annotations ${params.cohort}_annotations.ht \
            --output ${params.cohort}.mt

        # cut down on work folder space
        rm -r ${params.cohort}_annotations.ht

        tar --no-xattrs -czf ${params.cohort}.mt.tar.gz ${params.cohort}.mt
        rm -r ${params.cohort}.mt
        """
}
