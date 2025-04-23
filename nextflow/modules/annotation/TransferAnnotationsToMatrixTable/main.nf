
process TransferAnnotationsToMatrixTable {
    container params.container

    input:
        path(compressed_ht)
        tuple path(vcf), path(tbi)

    publishDir params.cohort_output_dir, mode: 'copy'

    output:
        path("${params.cohort}.mt.tar")

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
		tar --no-xattrs --remove-files -cf ${params.cohort}.mt.tar ${params.cohort}.mt
        """
}
