process MergeVcfsWithBcftools {
    container params.container

    input:
        tuple val(cohort), path(vcfs), path(tbis)
        path ref_genome

    output:
        tuple val(cohort), path("${cohort}_merged.vcf.bgz"), path("${cohort}_merged.vcf.bgz.tbi")

    script:
        def input = (vcfs.collect().size() > 1) ? vcfs.sort{ svcf -> svcf.name } : vcfs
        """
        set -euo pipefail

        # Each single-sample VCF is split and left-aligned before merging (populationgenomics/talos#729).
        # Merging first and splitting afterwards corrupts sites where an allele sits in one sample's
        # multi-allelic record but alone in another's (child C>CT,CTT 1/2, father C>CTT 0/1): the merge
        # pairs records by their full ALT list, so C>CTT comes out of the later split twice, each row
        # holding part of the genotypes - a false de novo on one, a lost allele on the other.
        # With every input biallelic and identically represented, the merge pairs each allele across samples.
        # -d exact drops records a caller wrote twice within one sample (seen in DRAGEN output)
        mkdir -p normalised
        for vcf in ${input.join(' ')}; do
            bcftools norm \
                -m -any \
                -f ${ref_genome} \
                -d exact \
                -Ob1 \
                --no-version \
                -o "normalised/\$(basename "\${vcf}").bcf" \
                -W \
                "\${vcf}"
        done

        # https://github.com/samtools/bcftools/issues/1189
        # "-m none" means don't merge multi-allelic sites, keep everything atomic
        # -0 to set all missing genotypes to HomWT - gap-filling with Missing (default) reduces the AN, so callset
        # frequency filters can appear to show an inflated AC/AN ratio
        bcftools merge \
        	--force-single \
        	-m none \
        	-0 \
        	-Oz \
        	--no-version \
        	-o "${cohort}_merged.vcf.bgz" \
        	-W=tbi \
        	normalised/*.bcf
        """
}
