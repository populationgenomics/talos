process AnnotatedVcfIntoMatrixTable {
    container params.container

    input:
        path vcf
        path gene_bed
        path mane

    publishDir params.cohort_output_dir

    output:
        path "${vcf.simpleName}_annotations.mt"

    script:
        """
        set -ex

        python -m talos.annotation_scripts.annotated_vcf_into_matrixtable \
            --input ${vcf} \
            --gene_bed ${gene_bed} \
            --output ${vcf.simpleName}_annotations.mt \
            --mane ${mane}

        # tidy up all checkpoints
        rm -r checkpoint*
        """
}
