
process ParseAlphaMissenseIntoHt {
    container params.container

    // parse AM data as a Hail Table
    publishDir params.processed_annotations, mode: 'copy'

    input:
        path am_tsv

    output:
        path "alphamissense_38.ht"

    script:
        """
        ParseAlphaMissenseIntoHt \
            --am_tsv ${am_tsv} \
            --ht_out alphamissense_38.ht
        rm temp.json
        """
}
