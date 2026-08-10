process ResummariseRawSubmissions {
    container params.container

    input:
        path variant_summary
        path submission_summary
        val timestamp

    output:
        tuple path("clinvarbitration_${timestamp}.vcf.bgz"), path("clinvarbitration_${timestamp}.vcf.bgz.tbi"), emit: "vcf"
        path "clinvarbitration_${timestamp}.ht", emit: "ht"
        path "clinvarbitration_${timestamp}.tsv", emit: "tsv"

    // Generates
    // clinvarbitration_XX.vcf.bgz + index - VCF containing only pathogenic SNV entries, feeds into annotation
    // clinvarbitration_XX.ht - a Hail Table containing the summarised data entries
    // clinvarbitration_XX.tsv - the same decisions as a TSV, consumed by the streaming (non-Hail) processes
    script:
        """
        set -euo pipefail

        python3 -m clinvarbitration.scripts.resummarise_clinvar \
            -v "${variant_summary}" \
            -s "${submission_summary}" \
            -o "clinvarbitration_${timestamp}" \
            -b "${params.clinvar_blacklist}"
        """
}
