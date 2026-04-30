process EncodeNapogee {
    container params.container

    input:
        path vcf

    output:
        path "napogee.zip"

    script:
        """
        set -euo pipefail

        echtvar encode napogee.zip /talos/echtvar/napogee_config.json ${vcf}
        """
}
