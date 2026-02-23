process MergeVcfsWithBcftools {
    container params.container

    publishDir "${params.outdir}/${cohort}_outputs"

    input:
        tuple val(cohort), path(vcfs), path(tbis)
        path ref_genome

    output:
        tuple val(cohort), path("${cohort}_merged.vcf.bgz")

    script:
        // should bump this check earlier tbh - if one VCF don't call this. Maybe just trust users
        def input = (vcfs.collect().size() > 1) ? vcfs.sort{ it.name } : vcfs
        """
        set -e
        # https://github.com/samtools/bcftools/issues/1189
        # "-m none" means don't merge multi-allelic sites, keep everything atomic, we're splitting in the next step
        # -0 to set all missing genotypes to HomWT - gap-filling with Missing (default) reduces the AN, so callset
        # frequency filters can appear to show an inflated AC/AN ratio
        bcftools merge \
        	--force-single \
        	-m none \
        	-0 \
        	-Ou \
        	--no-version \
        	-o "${params.cohort}_merged.vcf.bgz" \
        	$input
        """
}
