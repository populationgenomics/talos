process ReformatAnnotatedVcfIntoHailTable {
    container params.container

    input:
        path vcf
        path alphamissense
        path gene_bed
        path mane

    publishDir params.cohort_output_dir

    output:
        path "${params.cohort}_annotations.ht"

    script:
        """
        set -ex

        ReformatAnnotatedVcfIntoHailTable \
            --input ${vcf} \
            --am ${alphamissense} \
            --gene_bed ${gene_bed} \
            --output ${params.cohort}_annotations.ht \
            --mane ${mane}
        """
}
