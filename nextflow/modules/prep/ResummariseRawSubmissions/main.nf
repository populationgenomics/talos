process ResummariseRawSubmissions {
    container params.container

    input:
        path variant_summary
        path submission_summary
        val timestamp

    output:
        path "clinvarbitration_${timestamp}.vcf.bgz", emit: "vcf"
        path "clinvarbitration_${timestamp}.vcf.bgz.tbi", emit: "vcf_idx"
        path "clinvarbitration_${timestamp}.ht", emit: "ht"

    // Generates
    // clinvarbitration_XX.vcf.bgz + index - VCF containing only pathogenic SNV entries, feeds into annotation
    // clinvarbitration_XX.ht - a Hail Table containing the summarised data entries
    script:
    """
    python3 -m clinvarbitration.scripts.resummarise_clinvar \
        -v "${variant_summary}" \
        -s "${submission_summary}" \
        -o "clinvarbitration_${timestamp}" \
        -b "${params.clinvar_blacklist}"

    # remove the byproduct TSV which was read into a HT
    rm clinvarbitration_${timestamp}.tsv
    """
}
