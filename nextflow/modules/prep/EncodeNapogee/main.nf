process EncodeNapogee {
    container params.container
    publishDir params.processed_annotations, mode: 'copy'

    input:
        path tsv

    output:
        tuple path("napogee.vcf.gz"), path("napogee.zip")

    script:
        """
        set -euo pipefail

        python -m talos.scripts.parse_napogee \
            --input ${tsv} \
            --output napogee.vcf.gz

        echtvar encode napogee.zip /talos/echtvar/napogee_config.json napogee.vcf.gz
        """
}
