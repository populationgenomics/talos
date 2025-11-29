process EncodeAlphaMissense {
    container params.container

    publishDir params.processed_annotations, mode: 'copy'

    input:
        path vcf

    output:
        path "alphamissense.zip"

    script:
        """
        echtvar encode alphamissense.zip /talos/echtvar/am_config.json ${vcf}
        """
}
