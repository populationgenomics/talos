process ParseManeIntoJson {
    container params.container

    publishDir params.processed_annotations, mode: 'copy'

    input:
    	path mane_summary

    output:
        path "mane.json", emit: json

    script:
    """
    python -m talos.annotation_scripts.parse_mane_into_json \
        --input ${mane_summary} \
        --output mane.json \
        --format json
    """
}
