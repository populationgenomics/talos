process AnnotateSvWithSvafotate {
    container params.svafotate_container

    input:
        tuple val(cohort), path(vcf)
        path svafotate_bed

    output:
        tuple val(cohort), path("${cohort}_sv_popaf.vcf")

    script:
        // -s gnomAD is mandatory: without it Max_AF becomes the maximum across all four callsets in the BED,
        // which would be mislabelled as gnomAD. -f must be set explicitly - its 0.001 default matches any
        // single-basepair overlap. The output is uncompressed; this image carries no htslib CLI, so the
        // compression and indexing happen in RenameSvAfFields
        """
        set -euo pipefail

        svafotate annotate \
            -v ${vcf} \
            -o ${cohort}_sv_popaf.vcf \
            -b ${svafotate_bed} \
            -s gnomAD \
            -f ${params.sv_overlap_fraction} \
            -a best \
            --ins \
            --cpu ${task.cpus}
        """
}
