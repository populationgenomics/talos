
process UnifiedPanelAppParser {
    container params.container

    publishDir params.output_dir, mode: 'copy'

    input:
    	path talos_config
        path panelapp_cache
        path pedigree
        path hpo
        path check_file

    output:
        path "${params.cohort}_panelapp.json"

    """
    export TALOS_CONFIG=${talos_config}
    UnifiedPanelAppParser \
        --input $panelapp_cache \
        --output ${params.cohort}_panelapp.json \
        --pedigree $pedigree \
        --hpo $hpo
    """
}
