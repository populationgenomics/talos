process AnnotatedVcfIntoMatrixTable {
    container params.container

    input:
        tuple val(cohort), path(vcf)
        path gene_bed
        path mane

    publishDir "${params.outdir}/${cohort}_outputs", mode: 'copy'

    output:
        tuple val(cohort), path("${vcf.simpleName}_annotations.mt")

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
