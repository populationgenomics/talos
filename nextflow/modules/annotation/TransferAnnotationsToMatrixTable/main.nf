process TransferAnnotationsToMatrixTable {
    container params.container

    input:
        path annotations
        tuple path(vcf), path(tbi)

    publishDir params.cohort_output_dir, mode: 'copy'

    output:
        path "${params.cohort}.mt"

    script:
        """
        set -ex

        TransferAnnotationsToMatrixTable \
            --input ${vcf} \
            --annotations ${annotations} \
            --output ${params.cohort}.mt

        # cut down on work folder space
        rm -r ${params.cohort}_annotations.ht
        """
}
