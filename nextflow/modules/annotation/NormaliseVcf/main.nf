process NormaliseVcf {
    container params.container

    // take the merged VCF and index
    // ref_genome here is used to create parsimonious representations
    input:
        tuple val(cohort), path(vcf)
        path ref_genome

    output:
        tuple val(cohort), path("${vcf.simpleName}_merged_filtered.vcf.bgz"), path("${vcf.simpleName}_merged_filtered.vcf.bgz.tbi")

    script:
    """
    set -euo pipefail

    # multiple ways to reach this point, some do & some don't strictly index
    tabix ${vcf}
    bcftools norm \
    	-m -any \
    	-f ${ref_genome} \
        -Ou ${vcf} \
        --no-version | \
    bcftools +fill-tags \
        -Oz \
        --no-version \
        -o "${vcf.simpleName}_merged_filtered.vcf.bgz" \
        -W=tbi - -- -t AC,AF,AN
    """
}
