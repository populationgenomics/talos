
process HPOFlagging {
    container params.container

    // flag results in the result JSON which are good phenotypic matches
    publishDir params.output_dir, mode: 'copy'

    input:
        path talos_result_json
        path gene_symbol_map
        path gene_to_phenotype
        path phenio_db
        path talos_config

    output:
        path "${params.cohort}_pheno_annotated_report.json", emit: "pheno_annotated"
        path "${params.cohort}_pheno_filtered_report.json", emit: "pheno_filtered"

    // the command - first decompress the phenio db, then run the HPOFlagging script
    // then delete the decompressed db (it's huge)
    """
    gunzip -f -d ${phenio_db}
    export TALOS_CONFIG=${talos_config}
    HPOFlagging \
         --input ${talos_result_json} \
         --gene_map ${gene_symbol_map} \
         --gen2phen ${gene_to_phenotype} \
         --phenio phenio.db \
         --output ${params.cohort}_pheno_annotated_report.json \
         --phenout ${params.cohort}_pheno_filtered_report.json
    rm phenio.db
    """
}
