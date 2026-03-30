process HPOFlagging {
    container params.container

    input:
        tuple val(cohort), path(talos_result_json), path(talos_config)
        path gene_symbol_map
        path gene_to_phenotype
        path phenio_db
        val timestamp

    output:
        tuple val(cohort), path("${cohort}_full_report_${timestamp}.json")

    script:
        """
        set -euo pipefail

        export TALOS_CONFIG=${talos_config}
        python -m talos.hpo_flagging \
             --input ${talos_result_json} \
             --mane_json ${gene_symbol_map} \
             --gen2phen ${gene_to_phenotype} \
             --phenio ${phenio_db} \
             --output ${cohort}_full_report_${timestamp}.json
        """
}
