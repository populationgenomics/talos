process SortCpxIntervals {
    container params.container

    input:
        tuple val(cohort), path(vcf), path(tbi)

    // same VCF, with each complex variant's CPX_INTERVALS in coordinate order
    output:
        tuple val(cohort), path("${cohort}_cpx_sorted.vcf.bgz"), path("${cohort}_cpx_sorted.vcf.bgz.tbi")

    script:
        // GATK SVAnnotate assumes CPX_INTERVALS is coordinate-sorted and aborts the whole run on the first
        // record where it is not - GATK-SV writes delINVdel as DEL,DEL,INV, which trips it immediately.
        // Full reasoning lives in the script's docstring. This runs ahead of AnnotateSvWithGatk rather than
        // inside it because the sort needs the Talos container, not the GATK one.
        """
        set -euo pipefail

        python -m talos.annotation_scripts.sort_cpx_intervals \
            --input ${vcf} \
            --output ${cohort}_cpx_sorted.vcf

        bgzip -c ${cohort}_cpx_sorted.vcf > ${cohort}_cpx_sorted.vcf.bgz
        tabix -p vcf ${cohort}_cpx_sorted.vcf.bgz
        """
}
