process TransferAnnotationsToMatrixTable {
    container params.container

    input:
        path mt
        path annotations
        path regions

    publishDir params.cohort_output_dir, mode: 'copy'

    output:
        path "${params.cohort}.mt"

    script:
        """
        set -ex

        TransferAnnotationsToMatrixTable \
            --input ${mt} \
            --annotations ${annotations} \
            --output ${params.cohort}.mt \
            --regions ${regions}
        """
}
