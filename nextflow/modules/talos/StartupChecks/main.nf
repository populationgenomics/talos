process StartupChecks {
    container params.container

    input:
        tuple val(cohort), path(mts), path(pedigree), path(talos_config), path(history), path(ext), path(seqr)
        path clinvar

    output:
        tuple val(cohort), path("${cohort}_checked")

    script:
        def mt_string = (mts.collect().size() > 1) ? mts.sort{ it.name } : mts

        """
        set -e
        export TALOS_CONFIG=${talos_config}

        python -m talos.startup_checks \\
            --mt ${mt_string} \\
            --pedigree ${pedigree} \\
            --clinvar ${clinvar}

        echo "success" > "${cohort}_checked"
        """
}
