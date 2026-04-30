process EncodeAlphaMissense {
    container params.container

    input:
        path am_tsv

    output:
        path "alphamissense.zip"
        path "alphamissense.vcf.gz"

    script:
        """
        set -euo pipefail

        python -m talos.annotation_scripts.parse_alphamissense \
            --input ${am_tsv} \
            --output alphamissense.vcf.gz

        echtvar encode alphamissense.zip /talos/echtvar/am_config.json alphamissense.vcf.gz
        """
}
