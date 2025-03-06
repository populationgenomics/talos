
process localise_alphamissense {
    container params.hail_docker

    // parse AM data as a Hail Table
    publishDir params.generic_output_dir, mode: 'copy'

    output:
        path("alphamissense.tsv.gz")

    script:
        """
        wget ${params.alphamissense_url} -O alphamissense.tsv.gz
        """
}
