process EncodeMitimpact {
    container params.container

    input:
        path mitimpact_tsv

    output:
        path "mitimpact.vcf.gz"
        path "mitimpact.zip"

    script:
        """
        set -euo pipefail


        python -m talos.annotation_scripts.parse_mitimpact \
            --input ${am_tsv} \
            --output mitimpact.vcf.gz

        echtvar encode mitimpact.zip /talos/echtvar/mitimpact_config.json mitimpact.vcf.gz
        """
}
