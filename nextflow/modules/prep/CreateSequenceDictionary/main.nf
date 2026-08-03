process CreateSequenceDictionary {
    container params.gatk_container

    // This is the one process in the pipeline using publishDir rather than the workflow output block. The
    // dictionary belongs to the reference genome, not to a run's results, so params.ref_dict points into
    // processed_annotations - and the output block cannot write outside outputDir. Publishing it here means
    // the existence check in sv_annotation.nf skips this process on every subsequent run.
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
