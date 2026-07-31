process RenameSvAfFields {
    container params.container

    input:
        tuple val(cohort), path(vcf), path(talos_config)

    output:
        tuple val(cohort), path("${cohort}_sv_annotated.vcf.bgz"), path("${cohort}_sv_annotated.vcf.bgz.tbi")

    script:
        // the field prefix is read from the Talos config, the same key run_hail_filtering_sv uses, so the
        // written and expected field names cannot drift apart. SVAFotate wrote an uncompressed VCF, so bgzip
        // and index it here in the Talos container
        """
        set -euo pipefail

        export TALOS_CONFIG=${talos_config}

        python -m talos.annotation_scripts.rename_sv_af_fields \
            --input ${vcf} \
            --output ${cohort}_sv_renamed.vcf

        bgzip -c ${cohort}_sv_renamed.vcf > ${cohort}_sv_annotated.vcf.bgz
        tabix -p vcf ${cohort}_sv_annotated.vcf.bgz
        """
}
