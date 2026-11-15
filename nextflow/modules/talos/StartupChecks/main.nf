process StartupChecks {
    container params.container

    // check the config, pedigree, and the input contract of one annotated shard. Every shard is
    // written by the same processes over the same samples, so one is representative of all
    input:
        tuple val(cohort), path(vcf), path(pedigree), path(talos_config)

    output:
        tuple val(cohort), path("${cohort}_checked")

    script:
        """
        set -euo pipefail

        export TALOS_CONFIG=${talos_config}

        python -m talos.startup_checks \\
            --vcf ${vcf} \\
            --pedigree ${pedigree}

        echo "success" > "${cohort}_checked"
        """
}
