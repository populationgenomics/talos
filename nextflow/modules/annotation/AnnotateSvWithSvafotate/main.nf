process AnnotateSvWithSvafotate {
    container params.svafotate_container

    input:
        tuple val(cohort), path(vcf), path(tbi)
        path svafotate_bed

    // writes an uncompressed VCF - this image carries no htslib CLI, so compression and indexing
    // are left to RenameSvAfFields, which runs in the Talos container
    output:
        tuple val(cohort), path("${cohort}_sv_popaf.vcf")

    script:
        """
        set -euo pipefail

        # -s gnomAD is required. Without it Max_AF is the maximum across the BED's four merged callsets
        # (gnomAD, CCDG, TOPMed, ThousG), which would be mislabelled as a gnomAD frequency downstream.
        #
        # -f must be set explicitly. It defaults to 0.001, which is effectively "any single-basepair overlap
        # matches", and would annotate a 200bp deletion with a 3Mb gnomAD deletion's frequency.
        #
        # The reference BED is passed exactly as downloaded. SVAFotate strips the `chr` prefix from the query
        # VCF but not from this BED, so its Ensembl-style contig names are already correct - adding a `chr`
        # prefix here would silently zero every frequency in the callset.
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
