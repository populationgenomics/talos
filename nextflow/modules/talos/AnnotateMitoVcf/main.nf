process AnnotateMitoVcf {
    container params.container

    input:
        tuple val(cohort), path(vcf), path(panelapp_data), path(pedigree), path(talos_config)
        path ref_fa
        path gff3
        path clinvar
        path panelapp
        path mitimpact
        path mitotip
        path napogee

    output:
        tuple val(cohort), path("${cohort}_mito_labelled.vcf.bgz")

    script:
        """
        set -euo pipefail

        bcftools csq \
            --force \
            -f "${ref_fa}" \
            -g "${gff3}" \
            --local-csq \
            -C 2 \
            --threads 4 \
            -B 10 \
            --unify-chr-names 'chr,-,chr' \
            -Oz -o "${cohort}_mito_csq_annotated.vcf.bgz" \
            ${vcf}

        echtvar anno \
            -e ${napogee} \
            -e ${mitimpact} \
            -e ${mitotip} \
            ${cohort}_mito_csq_annotated.vcf.bgz" \
            ${cohort}_mito_all_annotated.vcf.bgz"

        export TALOS_CONFIG=${talos_config}

        python -m talos.reformat_and_label_mito_vcf \\
            --input "${cohort}_mito_all_annotated.vcf.bgz" \\
            --output "${cohort}_mito_labelled.vcf.bgz" \\
            --pedigree ${pedigree} \\
            --panelapp ${panelapp} \\
            --clinvar ${clinvar}
        """
}
