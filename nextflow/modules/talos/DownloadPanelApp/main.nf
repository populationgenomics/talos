def timestamp = new java.util.Date().format('yyyy-MM')

process DownloadPanelApp {
    container params.container

    publishDir params.processed_annotations, mode: 'copy'

    input:
        path mane
        path talos_config

    output:
        path "panelapp_${timestamp}.json", emit: json

    """
    export TALOS_CONFIG=${talos_config}
    DownloadPanelApp --output panelapp_${timestamp}.json --mane ${mane}
    """
}
