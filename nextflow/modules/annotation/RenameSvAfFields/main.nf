process RenameSvAfFields {
    container params.container

    input:
        tuple val(cohort), path(vcf), path(talos_config)

    // Max_AF -> {prefix}_sv_AF, Best_gnomAD_ID -> {prefix}_sv_SVID, prefix read from the Talos config
    output:
        tuple val(cohort), path("${cohort}_sv_annotated.vcf.bgz"), path("${cohort}_sv_annotated.vcf.bgz.tbi")

    script:
        """
        set -euo pipefail

        export TALOS_CONFIG=${talos_config}

        python -m talos.annotation_scripts.rename_sv_af_fields \
            --input ${vcf} \
            --output ${cohort}_sv_annotated.vcf

        bgzip -c ${cohort}_sv_annotated.vcf > ${cohort}_sv_annotated.vcf.bgz
        tabix -p vcf ${cohort}_sv_annotated.vcf.bgz
        """
}
