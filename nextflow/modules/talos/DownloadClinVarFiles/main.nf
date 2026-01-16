process DownloadClinVarFiles {
    container params.container

    publishDir params.large_files, mode: 'copy'

    def timestamp = new java.util.Date().format('yyyy-MM')

    output:
        path "submissions_${timestamp}.txt.gz", emit: submissions
        path "variants_${timestamp}.txt.gz", emit: variants

    """
    wget '${params.submission_summary}' -O submissions_${timestamp}.txt.gz &
    wget '${params.variant_summary}' -O variants_${timestamp}.txt.gz &
    wait
    """
}
