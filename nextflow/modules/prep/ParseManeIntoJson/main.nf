process ParseManeIntoJson {
    container params.container

    input:
    	path mane_summary

    output:
        path "mane.json"

    script:
    """
    ParseManeIntoJson --input ${mane_summary} --output mane.json --format json
    """
}
