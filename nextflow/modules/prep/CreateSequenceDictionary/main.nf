process CreateSequenceDictionary {
    container params.gatk_container

    publishDir params.processed_annotations, mode: 'copy'

    // large_files ships ref.fa and ref.fa.fai, but no sequence dictionary
    input:
        path ref_fa

    // GATK SVAnnotate requires a .dict for contig ordering
    output:
        path "ref.dict"

    script:
        """
        set -euo pipefail

        gatk CreateSequenceDictionary \
            -R ${ref_fa} \
            -O ref.dict
        """
}
