process RegisterAnalysis {
    container params.container
    tag "${cohort}"

    // record one published file as a Metamist analysis of `type`, linked to the Metamist cohort, in the cohort's own
    // project. `output` is the published URI, built by cpg_main.nf from the same analysis_dir publishDir copies to.
    // Included once per analysis type (RegisterResults, RegisterReport). Needs Metamist credentials in the task
    input:
        // analysis_meta: [key: value], each value without spaces
        tuple val(cohort), val(project), val(output), val(analysis_meta)
        val(type)
        // how long to wait for publishDir to finish copying the output
        val(wait_minutes)

    output:
        val(cohort)

    script:
        def attempts = wait_minutes * 6
        def meta_args = analysis_meta.collect { key, value -> "${key}=${value}" }.join(' ')
        """
        set -euo pipefail

        # publishDir copies asynchronously, so the output may not be in GCS yet
        for i in \$(seq ${attempts}); do
            python -c "import sys; from cpg_utils.existence_checks import exists; sys.exit(not exists(sys.argv[1]))" \\
                ${output} && break
            sleep 10
        done

        # fails if the output is still missing, or isn't in GCS
        python -m cpg_utils.metamist_registration \\
            --project ${project} \\
            --output ${output} \\
            --type ${type} \\
            --cohorts ${cohort} \\
            --meta ${meta_args}
        """
}
