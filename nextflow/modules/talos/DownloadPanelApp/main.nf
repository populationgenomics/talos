
process DownloadPanelApp {
    container params.container

    publishDir params.processed_annotations, mode: 'copy'

    input:
        path mane
        path talos_config

    output:
        path "panelapp_download.json"

    """
    export TALOS_CONFIG=${talos_config}
    DownloadPanelApp --output panelapp_download.json --mane ${mane}
    """
}
