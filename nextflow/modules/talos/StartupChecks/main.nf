process StartupChecks {
    container params.container

    publishDir params.output_dir, mode: 'copy'

    input:
        path mts
        path pedigree
        path clinvar
        path talos_config

    output:
        path "${params.cohort}_checked"

    script:
        def mt_string = (mts.collect().size() > 1) ? mts.sort{ it.name } : mts

        """
        set -e
        export TALOS_CONFIG=${talos_config}

        python -m talos.startup_checks \\
            --mt ${mt_string} \\
            --pedigree ${pedigree} \\
            --clinvar ${clinvar}

        echo "success" > "${params.cohort}_checked"
        """
}
