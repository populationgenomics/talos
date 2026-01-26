process DownloadClinVarFiles {
    container params.container

    input:
        val timestamp

    output:
        path "submissions_${timestamp}.txt.gz", emit: submissions
        path "variants_${timestamp}.txt.gz", emit: variants

    script:
    """
    wget '${params.submission_summary}' -O submissions_${timestamp}.txt.gz &
    wget '${params.variant_summary}' -O variants_${timestamp}.txt.gz &
    wait
    """
}