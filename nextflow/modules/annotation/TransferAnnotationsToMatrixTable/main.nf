
process TransferAnnotationsToMatrixTable {
    container params.container

    input:
        path annotations
        tuple path(vcf), path(tbi)

    publishDir params.cohort_output_dir, mode: 'move'

    output:
        path "${params.cohort}.mt"

    script:
        """
        set -ex

        TransferAnnotationsToMatrixTable \
            --input ${vcf} \
            --annotations ${annotations} \
            --output ${params.cohort}.mt
        """
}
