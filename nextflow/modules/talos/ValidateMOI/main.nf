
process ValidateMOI {
    container params.container


    // process the labelled variants
    publishDir params.output_dir, mode: 'copy'

    input:
        tuple path(labelled_vcf), path(labelled_vcf_index)
        path panelapp
        path pedigree
        path talos_config
        path history_db

	def timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm')

    output:
        path "${params.cohort}_results.json", emit: json
        path "${timestamp}.db", emit: db

    """
    export TALOS_CONFIG=${talos_config}

    # if the db file exists, copy it
    if [ -f "${history_db}" ]; then
        cp "${history_db}" "talos_history_${timestamp}.db"
    fi

    ValidateMOI \
        --labelled_vcf ${labelled_vcf} \
        --panelapp ${panelapp} \
        --pedigree ${pedigree} \
        --output ${params.cohort}_results.json \
        --db "talos_history_${timestamp}.db"
    """
}
