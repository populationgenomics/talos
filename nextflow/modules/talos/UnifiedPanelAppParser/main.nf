
process UnifiedPanelAppParser {
    container params.container

    publishDir "${params.outdir}/${cohort}_outputs", mode: 'copy'

    input:
        tuple val(cohort), path(check_file)
    	path talos_config
        path panelapp_cache
        path pedigree
        path hpo

    output:
        tuple val(cohort), path("${cohort}_panelapp.json")

    """
    export TALOS_CONFIG=${talos_config}
    python -m talos.unified_panelapp_parser \
        --input $panelapp_cache \
        --output ${cohort}_panelapp.json \
        --pedigree $pedigree \
        --hpo $hpo
    """
}
