
process FilterVcfToBedWithBcftools {
    container params.container

    // take the merged VCF and index
    // also take a BED file of regions to focus analysis on/filter VCF to
    // ref_genome here is used to create parsimonious representations
    input:
        path vcf
        path tbi
        path bed_file
        path ref_genome

    publishDir params.cohort_output_dir, mode: 'copy'

    output:
        tuple \
            path("${params.cohort}_merged_filtered.vcf.bgz"), \
            path("${params.cohort}_merged_filtered.vcf.bgz.tbi")

    script:
    """
    bcftools norm \
    	-m -any \
    	-f ${ref_genome} \
        --write-index=tbi \
        -R ${bed_file} \
        -Oz -o "${params.cohort}_merged_filtered.vcf.bgz" \
        ${vcf}
    """
}
