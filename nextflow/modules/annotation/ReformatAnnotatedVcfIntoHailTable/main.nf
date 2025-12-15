process ReformatAnnotatedVcfIntoHailTable {
    container params.container

    input:
        path vcf
        path gene_bed
        path mane

    publishDir params.cohort_output_dir

    output:
        path "${vcf.baseName}_annotations.ht"

    script:
        """
        set -ex

        ReformatAnnotatedVcfIntoHailTable \
            --input ${vcf} \
            --gene_bed ${gene_bed} \
            --output ${vcf.baseName}_annotations.ht \
            --mane ${mane}
        """
}
