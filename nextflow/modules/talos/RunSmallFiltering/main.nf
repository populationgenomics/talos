process RunSmallFiltering {
    container params.container

    // filter and category-label one annotated shard in a single cyvcf2 pass. The MANE JSON is
    // consumed here rather than in the annotation workflow, because the consequence reshaping and
    // the labelling are now the same pass
    input:
        tuple val(cohort), path(vcf), path(panelapp_data), path(pedigree), path(talos_config)
        path mane
        path clinvar_pm5

    output:
        tuple val(cohort), path("${vcf.simpleName}_labelled.vcf.bgz"), path("${vcf.simpleName}_labelled.vcf.bgz.tbi")

    script:
        """
        set -euo pipefail

        export TALOS_CONFIG=${talos_config}

        python -m talos.run_small_filtering \
            --input ${vcf} \
            --panelapp ${panelapp_data} \
            --pedigree ${pedigree} \
            --mane ${mane} \
            --pm5 ${clinvar_pm5} \
            --output "${vcf.simpleName}_labelled.vcf.bgz"
        tabix "${vcf.simpleName}_labelled.vcf.bgz"
        """
}
