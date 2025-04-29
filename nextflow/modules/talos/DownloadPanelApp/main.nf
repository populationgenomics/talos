
process DownloadPanelApp {
    container params.container

    publishDir params.processed_annotations, mode: 'copy'

    input:
        path mane

    output:
        path "panelapp_download.json"

    """
    DownloadPanelApp --output panelapp_download.json --mane ${mane}
    """
}
