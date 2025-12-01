process AnnotateWithEchtvar {
    container params.container

    input:
        path vcf
        path gnomad_zip
        path am_zip

    // annotate this VCF with gnomAD & AlphaMissense data
    publishDir params.cohort_output_dir

    output:
        path("${vcf.baseName}_echtvar.vcf.bgz")

    // annotate with both sources, hard filter to 5% gnomAD 4 joint-exome/genome AF
    script:
        """
        set -ex
        echtvar anno \
            -e ${gnomad_zip} \
            -e ${am_zip} \
            -i "gnomad_AF_joint < 0.05" \
            ${vcf} \
            "${vcf.baseName}_echtvar.vcf.bgz"
        """
}
