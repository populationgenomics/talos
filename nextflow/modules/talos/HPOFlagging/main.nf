process HPOFlagging {
    container params.container

    // flag results in the result JSON which are good phenotypic matches
    publishDir "${params.outdir}/${cohort}_outputs", mode: 'copy'

    input:
        tuple val(cohort), path(talos_result_json), path(talos_config)
        path gene_symbol_map
        path gene_to_phenotype
        path phenio_db

	def timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm')

    output:
        tuple val(cohort), path("${cohort}_full_report_${timestamp}.json")

    """
    export TALOS_CONFIG=${talos_config}
    python -m talos.hpo_flagging \
         --input ${talos_result_json} \
         --mane_json ${gene_symbol_map} \
         --gen2phen ${gene_to_phenotype} \
         --phenio ${phenio_db} \
         --output ${cohort}_full_report_${timestamp}.json
    """
}
