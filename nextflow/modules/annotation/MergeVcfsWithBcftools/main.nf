
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
        # "-m none" means don't merge multi-allelic sites, keep everything atomic, we're splitting in the next step
        # -0 to set all missing genotypes to HomWT - gap-filling with Missing (default) reduces the AN, so callset
        # frequency filters can appear to show an inflated AC/AN ratio
        bcftools merge \
        	--force-single \
        	-m none \
        	-R ${regions} \
        	-0 \
        	-Ou \
        	$input | \
        bcftools norm \
			-m -any \
			-Ov \
			-f ${ref_genome} \
			-o temp.vcf \
			-
        bcftools +fill-tags temp.vcf -Oz -o "${params.cohort}_merged.vcf.bgz" -W=tbi -- -t AF
        rm temp.vcf
        """
}
