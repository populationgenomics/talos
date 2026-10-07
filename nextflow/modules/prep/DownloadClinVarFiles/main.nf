process DownloadClinVarFiles {
    container params.container

    errorStrategy {'retry'}
    maxRetries 3

    input:
        val timestamp
        val evidence_date

    output:
        path "submissions_${timestamp}.txt.gz", emit: submissions
        path "variants_${timestamp}.txt.gz", emit: variants

    script:
        def historical = evidence_date ? true : false
        def submission_url = historical ? "${params.clinvar_archive}/submission_summary_${timestamp}.txt.gz" : params.submission_summary
        def variant_url = historical ? "${params.clinvar_archive}/variant_summary_${timestamp}.txt.gz" : params.variant_summary
        def validation
        if (historical) {
            validation = """
                gzip -t submissions_${timestamp}.txt.gz
                gzip -t variants_${timestamp}.txt.gz
                sha256sum submissions_${timestamp}.txt.gz variants_${timestamp}.txt.gz
            """
        } else {
            validation = """
                wget '${submission_url}.md5' -O submissions.md5
                wget '${variant_url}.md5' -O variants.md5
                sed -E 's/[[:space:]].*/  submissions_${timestamp}.txt.gz/' submissions.md5 > check.md5
                sed -E 's/[[:space:]].*/  variants_${timestamp}.txt.gz/' variants.md5 >> check.md5
                md5sum --check check.md5
            """
        }
        """
        set -euo pipefail

        wget '${submission_url}' -O submissions_${timestamp}.txt.gz
        wget '${variant_url}' -O variants_${timestamp}.txt.gz
        ${validation}
        """
}
