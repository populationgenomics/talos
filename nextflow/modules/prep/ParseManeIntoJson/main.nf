process ParseManeIntoJson {
    container params.container

    input:
    	path mane_summary

    output:
        path "mane.json"

    script:
    """
    python -m talos.annotation_scripts.parse_mane_into_json \
        --input ${mane_summary} \
        --output mane.json \
        --format json
    """
}
