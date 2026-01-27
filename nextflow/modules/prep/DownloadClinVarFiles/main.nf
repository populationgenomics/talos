process DownloadClinVarFiles {
    container params.container

    input:
        val timestamp

    output:
        path "submissions_${timestamp}.txt.gz", emit: submissions
        path "variants_${timestamp}.txt.gz", emit: variants

    shell:
    """
    wget '!{params.submission_summary}' -O submissions_!{timestamp}.txt.gz
    wget '!{params.variant_summary}' -O variants_!{timestamp}.txt.gz
    wget '!{params.submission_summary}.md5' -O submissions.md5
    wget '!{params.variant_summary}'.md5 -O variants.md5
    awk '{print \$1, "submissions_!{timestamp}.txt.gz"}' submissions.md5 > check.md5
    awk '{print \$1, "variants_!{timestamp}.txt.gz"}' variants.md5 >> check.md5
    md5sum --check check.md5
    """
}
