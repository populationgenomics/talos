
process generate_roi {
    container params.hail_docker

    // generate the ROI file
    publishDir params.generic_output_dir, mode: 'copy'

    output:
        path "GRCh38.bed", emit: bed
        path "GRCh38.gff3.gz", emit: gff3

    script:
        def url_template = "${params.ensembl_gtf}".replaceAll("VER", "${params.ensembl_version}")
        """
        wget ${url_template} -O GRCh38.gff3.gz
        python3 /talos/generate_gene_roi.py \
            --gff3 GRCh38.gff3.gz \
            --output GRCh38.bed
        """
}
