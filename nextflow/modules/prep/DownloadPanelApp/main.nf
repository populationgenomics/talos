process DownloadPanelApp {
    container params.container

    input:
        path mane
        val timestamp

    output:
        path "panelapp_${timestamp}.json"

    script:
    """
    DownloadPanelApp --output panelapp_${timestamp}.json --mane ${mane}
    """
}

