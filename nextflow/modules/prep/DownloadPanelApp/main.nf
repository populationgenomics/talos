process DownloadPanelApp {
    container params.container
    publishDir params.processed_annotations, mode: 'copy'

    input:
        val timestamp

    output:
        path "panelapp_${timestamp}.json"

    script:
        """
        set -euo pipefail

        python -m talos.download_panelapp --output panelapp_${timestamp}.json
        """
}
