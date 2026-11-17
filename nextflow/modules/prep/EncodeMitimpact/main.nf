process EncodeMitimpact {
    container params.container
    publishDir params.processed_annotations, mode: 'copy'

    input:
        path tsv

    output:
        tuple path("mitimpact.vcf.gz"), path("mitimpact.zip")

    script:
        """
        set -euo pipefail

        python -m talos.scripts.parse_mitimpact \
            --input ${tsv} \
            --output mitimpact.vcf.gz

        echtvar encode mitimpact.zip /talos/echtvar/mitimpact_config.json mitimpact.vcf.gz
        """
}
