process ParseAlphaMissense {
    container params.container

    // parse AM data as a VCF, then encode as a condensed zip file
    publishDir params.processed_annotations, mode: 'copy'

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
