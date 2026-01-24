process NormaliseAndRegionFilterVcf {
    container params.container

    // take the merged VCF and index
    // also take a BED file of regions to focus analysis on/filter VCF to
    // ref_genome here is used to create parsimonious representations
    input:
        path vcf
        path bed_file
        path ref_genome

    publishDir params.cohort_output_dir, mode: 'copy'

    output:
        tuple \
            path("${vcf.simpleName}_merged_filtered.vcf.bgz"), \
            path("${vcf.simpleName}_merged_filtered.vcf.bgz.tbi")

    script:
    """
    tabix ${vcf}
    bcftools norm \
    	-m -any \
    	-f ${ref_genome} \
        -R ${bed_file} \
        -Ou ${vcf} \
        --no-version | \
    bcftools +fill-tags \
        -Oz \
        --no-version \
        -o "${vcf.simpleName}_merged_filtered.vcf.bgz" \
        -W=tbi - -- -t AC,AF,AN
    """
}
