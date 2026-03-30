process EncodeAlphaMissense {
    container params.container

    input:
        path vcf

    output:
        path "alphamissense.zip"

    script:
        """
        set -euo pipefail

        echtvar encode alphamissense.zip /talos/echtvar/am_config.json ${vcf}
        """
}
