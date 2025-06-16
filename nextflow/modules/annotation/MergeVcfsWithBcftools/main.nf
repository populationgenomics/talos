
process MergeVcfsWithBcftools {
    container params.container

    input:
        path vcfs
        path tbis
        path regions
        path ref_genome

    // merge the VCFs into a single VCF
    publishDir params.cohort_output_dir, mode: 'copy'

    output:
        tuple \
            path("${params.cohort}_merged.vcf.bgz"), \
            path("${params.cohort}_merged.vcf.bgz.tbi")

    script:

        def input = (vcfs.collect().size() > 1) ? vcfs.sort{ it.name } : vcfs
        """
        # https://github.com/samtools/bcftools/issues/1189
        # -m none means don't merge multi-allelic sites, keep everything atomic
        bcftools merge \
        	--force-single \
        	-m none \
        	-R ${regions} \
        	$input | \
        bcftools norm \
			-m -any \
			-f ${ref_genome} | \
        bcftools +fill-tags -Oz -o "${params.cohort}_merged.vcf.bgz" --write-index=tbi -- -t AF
        """
}
