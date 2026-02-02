process ParseAlphaMissense {
    container params.container

    input:
        path am_tsv

    output:
        path "alphamissense.vcf.gz"

    script:
        """
        python -m talos.annotation_scripts.parse_alphamissense \
            --input ${am_tsv} \
            --output alphamissense.vcf.gz
        """
}
