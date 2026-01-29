process DownloadPanelApp {
    container params.container

    input:
        path mane
        val timestamp

    output:
        path "panelapp_${timestamp}.json"

    script:
    """
    python -m talos.download_panelapp \
        --output panelapp_${timestamp}.json \
        --mane ${mane}
    """
}
