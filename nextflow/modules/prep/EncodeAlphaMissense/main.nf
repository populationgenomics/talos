process EncodeAlphaMissense {
    container params.container

    input:
        path tsv

    output:
        tuple path("alphamissense.vcf.gz"), path("alphamissense.zip")

    script:
        """
        set -euo pipefail

        python -m talos.annotation_scripts.parse_alphamissense \
            --input ${tsv} \
            --output alphamissense.vcf.gz

        echtvar encode alphamissense.zip /talos/echtvar/am_config.json alphamissense.vcf.gz
        """
}
