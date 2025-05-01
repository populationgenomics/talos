
process TransferAnnotationsToMatrixTable {
    container params.container

    input:
        path(compressed_ht)
        tuple path(vcf), path(tbi)

    publishDir params.cohort_output_dir, mode: 'link'

    output:
        path("${params.cohort}.mt")

    script:
        """
        set -ex
        tar --no-same-owner -xf ${compressed_ht}
        TransferAnnotationsToMatrixTable \
            --input ${vcf} \
            --annotations ${params.cohort}_annotations.ht \
            --output ${params.cohort}.mt

        # cut down on work folder space
        rm -r ${params.cohort}_annotations.ht
        """
}
