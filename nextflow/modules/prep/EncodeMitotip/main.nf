process EncodeMitotip {
    container params.container
    publishDir params.processed_annotations, mode: 'copy'

    input:
        path tsv

    output:
        tuple path("mitotip.vcf.gz"), path("mitotip.zip")

    script:
        """
        set -euo pipefail

        python -m talos.scripts.parse_mitotip \
            --input ${tsv} \
            --output mitotip.vcf.gz

        echtvar encode mitotip.zip /talos/echtvar/mitotip_config.json mitotip.vcf.gz
        """
}
