process UnifiedPanelAppParser {
    container params.container

    input:
        tuple val(cohort), path(check_file), path(talos_config), path(pedigree)
        path panelapp_cache
        path hpo
        val evidence_date

    output:
        tuple val(cohort), path("${cohort}_panelapp.json")

    script:
        def evidence_date_arg = evidence_date ? "--evidence-date ${evidence_date}" : ""
        """
        set -euo pipefail

        export TALOS_CONFIG=${talos_config}
        python -m talos.unified_panelapp_parser \
            --input $panelapp_cache \
            --output ${cohort}_panelapp.json \
            --pedigree $pedigree \
            --hpo $hpo ${evidence_date_arg}
        """
}
