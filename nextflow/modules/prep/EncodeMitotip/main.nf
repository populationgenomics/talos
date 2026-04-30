process EncodeMitotip {
    container params.container

    input:
        path vcf

    output:
        path "mitotip.zip"

    script:
        """
        set -euo pipefail

        echtvar encode mitotip.zip /talos/echtvar/mitotip_config.json ${vcf}
        """
}
