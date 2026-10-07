process DownloadPanelApp {
    container params.container

    input:
        val panelapp_key
        val evidence_date

    output:
        path "panelapp_${panelapp_key}.json"

    script:
        def evidence_date_arg = evidence_date ? "--evidence-date ${evidence_date}" : ""
        """
        set -euo pipefail

        python -m talos.download_panelapp \
            --output panelapp_${panelapp_key}.json ${evidence_date_arg}
        """
}
