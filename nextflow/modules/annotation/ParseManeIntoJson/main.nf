
process ParseManeIntoJson {
    container params.container

    publishDir params.processed_annotations, mode: 'copy'

    input:
    	path mane_summary

    output:
        path "mane.json", emit: json

    script:
    """
    ParseManeIntoJson --input ${mane_summary} --output mane.json --format json
    """
}
