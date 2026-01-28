process StartupChecks {
    container params.container

    publishDir params.output_dir, mode: 'copy'

    input:
        path mt
        path pedigree
        path clinvar
        path talos_config

    output:
        path "${params.cohort}_checked"

    // untar the ClinvArbitration data directory so we can check its contents
    """
    set -e
    export TALOS_CONFIG=${talos_config}

	python -m talos.startup_checks \\
        --mt ${mt} \\
        --pedigree ${pedigree} \\
        --clinvar ${clinvar}

    echo "success" > "${params.cohort}_checked"
    """
}
