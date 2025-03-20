
process ReformatAnnotatedVcfIntoHailTable {
    container params.container

    input:
        path(vcf)
        path(alphamissense)
        path(gene_bed)
        path(mane)

    publishDir params.cohort_output_dir, mode: 'copy'

    output:
        path("${params.cohort}_annotations.ht.tar.gz")

    // TODO write this script...
    script:
        """
        set -ex
        tar -xf ${alphamissense}
        ReformatAnnotatedVcfIntoHailTable \
            --input ${vcf} \
            --am alphamissense_38.ht \
            --gene_bed ${gene_bed} \
            --output ${params.cohort}_annotations.ht \
            --mane ${mane}

        # cut down on work folder space
        rm -r alphamissense_38.ht

        tar --no-xattrs -cf ${params.cohort}_annotations.ht.tar ${params.cohort}_annotations.ht

        rm -r ${params.cohort}_annotations.ht
        """
}
