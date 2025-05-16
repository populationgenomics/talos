
process ReformatAnnotatedVcfIntoHailTable {
    container params.container

    input:
        path vcf
        path alphamissense
        path gene_bed
        path mane

    publishDir params.cohort_output_dir, mode: 'copy'

    output:
        path "${params.cohort}_annotations.ht.tar"

    script:
        """
        set -ex
        tar --no-same-owner -xf ${alphamissense}
        ReformatAnnotatedVcfIntoHailTable \
            --input ${vcf} \
            --am alphamissense_38.ht \
            --gene_bed ${gene_bed} \
            --output ${params.cohort}_annotations.ht \
            --mane ${mane}

		tar --remove-files -cf ${params.cohort}_annotations.ht.tar ${params.cohort}_annotations.ht

        # cut down on work folder space
        rm -r alphamissense_38.ht
        """
}
