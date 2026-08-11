process ResummariseRawSubmissions {
    container params.container

    input:
        path variant_summary
        path submission_summary
        val timestamp

    output:
        tuple path("clinvarbitration_${timestamp}.vcf.bgz"), path("clinvarbitration_${timestamp}.vcf.bgz.tbi"), emit: "vcf"
        tuple path("clinvarbitration_${timestamp}.all.vcf.bgz"), path("clinvarbitration_${timestamp}.all.vcf.bgz.tbi"), emit: "all_vcf"
        path "clinvarbitration_${timestamp}.tsv", emit: "tsv"

    // Generates
    // clinvarbitration_XX.vcf.bgz + index - VCF containing only pathogenic SNV entries, feeds into the PM5 codon map
    // clinvarbitration_XX.all.vcf.bgz + index - every decision, Benign included. Encoded for echtvar at run time
    // clinvarbitration_XX.tsv - the same decisions as a TSV, consumed by the mito labelling process
    script:
        """
        set -euo pipefail

        python3 -m clinvarbitration.scripts.resummarise_clinvar \
            -v "${variant_summary}" \
            -s "${submission_summary}" \
            -o "clinvarbitration_${timestamp}" \
            -b "${params.clinvar_blacklist}" \
            --all_vcf "clinvarbitration_${timestamp}.all.vcf.bgz"
        """
}
