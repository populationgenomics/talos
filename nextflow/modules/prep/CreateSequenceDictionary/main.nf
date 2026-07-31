process CreateSequenceDictionary {
    container params.gatk_container

    // land the dictionary next to the reference in processed_annotations, so the sv_annotation.nf gate
    // (which checks params.ref_dict) finds it on a subsequent run and skips this pass over the 3GB FASTA
    publishDir params.processed_annotations, mode: 'copy', overwrite: true

    input:
        path ref_fa

    output:
        path 'ref.dict'

    script:
        """
        set -euo pipefail

        gatk CreateSequenceDictionary \
            -R ${ref_fa} \
            -O ref.dict
        """
}
