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

    tar --no-same-owner -zxf ${clinvar}

	python -m talos.StartupChecks \\
        --mt ${mt} \\
        --pedigree ${pedigree} \\
        --clinvar clinvarbitration_data/clinvar_decisions.ht \\
                  clinvarbitration_data/clinvar_decisions.pm5.ht

    echo "success" > "${params.cohort}_checked"
    rm -r clinvarbitration_data
    """
}
