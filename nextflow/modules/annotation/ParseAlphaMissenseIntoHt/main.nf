
process ParseAlphaMissenseIntoHt {
    container params.container

    // parse AM data as a Hail Table
    publishDir params.generic_output_dir, mode: 'copy'

    input:
        path(am_tsv)
        path(mane_json)

    output:
        path("alphamissense_isoforms_38.ht.tar")

    script:
        """
        ParseAlphaMissenseIntoHt \
            --am_tsv ${am_tsv} \
            --ht_out alphamissense_isoforms_38.ht \
            --mane_json ${mane_json}

        tar --no-xattrs -cf alphamissense_isoforms_38.ht.tar alphamissense_isoforms_38.ht
        rm -r alphamissense_isoforms_38.ht
        """
}
