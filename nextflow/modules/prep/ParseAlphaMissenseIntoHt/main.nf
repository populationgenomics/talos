process ParseAlphaMissenseIntoHt {
    container params.container

    input:
        path am_tsv

    output:
        path "alphamissense_38.ht"

    script:
        """
        python -m talos.annotation_scripts.parse_alphamissense \
            --am_tsv ${am_tsv} \
            --ht_out alphamissense_38.ht
        rm temp.json
        """
}
