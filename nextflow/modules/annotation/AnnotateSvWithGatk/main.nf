process AnnotateSvWithGatk {
    container params.gatk_container

    input:
        tuple val(cohort), path(vcf), path(vcf_tbi)
        path mane_gtf
        path noncoding_bed
        path ref_dict

    output:
        tuple val(cohort), path("${cohort}_sv_gatk.vcf.gz")

    script:
        """
        set -euo pipefail

        # SVAnnotate has no codec for a gzipped GTF, and fails with "no suitable codecs found", so decompress
        # it to a scratch file first - uncompressed, it is read with EnsemblGtfCodec
        gzip -dc ${mane_gtf} > mane.gtf

        gatk SVAnnotate \
            -V ${vcf} \
            --protein-coding-gtf mane.gtf \
            --non-coding-bed ${noncoding_bed} \
            --sequence-dictionary ${ref_dict} \
            -O ${cohort}_sv_gatk.vcf.gz
        """
}
