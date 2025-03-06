
process pull_and_parse_mane {
    container params.hail_docker

    publishDir params.generic_output_dir, mode: 'copy'

    output:
        path('mane.json')

    script:
    """
    wget ${params.mane} -O mane_summary.txt.gz
    python3 /talos/reformat_mane_summary.py --input mane_summary.txt.gz --output mane.json --format json
    """
}
