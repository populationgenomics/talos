process StartupChecks {
    container params.container

    input:
        tuple val(cohort), path(mts), path(pedigree), path(talos_config), path(history), path(ext), path(seqr), path(mito)
        path clinvar
        val evidence_date

    output:
        tuple val(cohort), path("${cohort}_checked")

    script:
        def mt_string = (mts.collect().size() > 1) ? mts.sort{ it.name } : mts
        def evidence_date_arg = evidence_date ? "--evidence-date ${evidence_date}" : ""

        """
        set -euo pipefail

        export TALOS_CONFIG=${talos_config}

        python -m talos.startup_checks \\
            --mt ${mt_string} \\
            --pedigree ${pedigree} \\
            --clinvar ${clinvar} ${evidence_date_arg}

        echo "success" > "${cohort}_checked"
        """
}
