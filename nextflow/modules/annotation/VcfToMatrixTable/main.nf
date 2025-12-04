process VcfToMatrixTable {
    container params.container

    input:
        path vcf
        path bed

    publishDir params.cohort_output_dir, mode: 'copy'

    output:
        path "${params.cohort}_initial.mt"

    script:
        """
        python -m talos.annotation_scripts.VcfToMatrixTable \
            --input ${vcf} \
            --regions ${bed} \
            --output ${params.cohort}_initial.mt
        rm -r ${params.cohort}_initial.mt._temp
        """
}
