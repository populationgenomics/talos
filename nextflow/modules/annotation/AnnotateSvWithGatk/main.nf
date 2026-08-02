process AnnotateSvWithGatk {
    container params.gatk_container

    input:
        tuple val(cohort), path(vcf), path(tbi)
        path mane_gtf
        path noncoding_bed
        path ref_dict

    // adds the PREDICTED_* consequence annotations, PREDICTED_LOF being the one Talos requires
    output:
        tuple val(cohort), path("${cohort}_sv_consequences.vcf.gz"), path("${cohort}_sv_consequences.vcf.gz.tbi")

    script:
        """
        set -euo pipefail

        # SVAnnotate has no codec for a gzipped GTF, so it has to be decompressed first.
        # No further preparation is needed - the MANE GTF is one MANE Select transcript per gene, which
        # satisfies SVAnnotate's canonical-transcript constraint as downloaded.
        gzip -dc ${mane_gtf} > mane.gtf

        # --non-coding-bed adds PREDICTED_NONCODING_BREAKPOINT and PREDICTED_NONCODING_SPAN. Talos reads
        # neither today, so these are carried for provenance and future use only - PREDICTED_LOF is unaffected.
        gatk SVAnnotate \
            -V ${vcf} \
            -O ${cohort}_sv_consequences.vcf.gz \
            --protein-coding-gtf mane.gtf \
            --non-coding-bed ${noncoding_bed} \
            --sequence-dictionary ${ref_dict}
        """
}
