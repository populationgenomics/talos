
process vcf_to_mt {
    container params.hail_docker

    input:
        path(vcf)
        path(alphamissense)
        path(gene_bed)
        path(mane)

    publishDir params.cohort_output_dir, mode: 'copy'

    output:
        path("${params.cohort}.mt.tar.gz")

    script:
        """
        set -ex
        tar -xf ${alphamissense}
        python3 /talos/convert_vcf_to_mt.py \
            --input ${vcf} \
            --am alphamissense.ht \
            --gene_bed ${gene_bed} \
            --output ${params.cohort}.mt \
            --mane ${mane}
        # cut down on work folder space
        rm -r alphamissense.ht

        tar -czf ${params.cohort}.mt.tar.gz ${params.cohort}.mt
        # not deleting the MT at the moment, someone might want to use it without decompressing the process output
        """
}
