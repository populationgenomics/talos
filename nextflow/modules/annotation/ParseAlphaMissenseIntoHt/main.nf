
process ParseAlphaMissenseIntoHt {
    container params.container

    // parse AM data as a Hail Table
    publishDir params.generic_output_dir, mode: 'copy'

    input:
        path(am_tsv)

    output:
        path("alphamissense.ht.tar.gz")

    script:
        """
        ParseAlphaMissenseIntoHt \
            --am_tsv ${am_tsv} \
            --ht_out alphamissense.ht
        tar -czf alphamissense.ht.tar.gz alphamissense.ht
        rm -r alphamissense.ht
        """
}
