process HPOFlagging {
    container params.container

    input:
        tuple val(cohort), path(talos_result_json), path(panelapp), path(talos_config)
        path gene_to_phenotype
        path phenio_db
        val timestamp

    output:
        tuple val(cohort), path("${cohort}_full_results_${timestamp}.json")

    script:
        """
        set -euo pipefail

        export TALOS_CONFIG=${talos_config}
        python -m talos.hpo_flagging \
             --input ${talos_result_json} \
             --panelapp ${panelapp} \
             --gen2phen ${gene_to_phenotype} \
             --phenio ${phenio_db} \
             --output ${cohort}_full_results_${timestamp}.json
        """
}
