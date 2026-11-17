process InputTsvFromMetamist {
    container params.container

    // query Metamist for each cohort's inputs, generate its pedigree, external IDs, seqr map and config, and write
    // those straight to <cohort root>/analysis_<run_date>/inputs/ - so nothing is published from here. The TSV
    // emitted holds full URIs to all of them, one row per cohort. Needs Metamist credentials in the task environment
    input:
        val(cohort_ids)
        val(access_level)
        val(sequencing_type)
        path(talos_config)
        val(run_date)

    output:
        path("input.tsv")

    script:
        """
        set -euo pipefail

        python -m talos.cpg_internal.input_tsv_from_metamist \\
            --cohort_ids ${cohort_ids.join(' ')} \\
            --access_level ${access_level} \\
            --sequencing_type ${sequencing_type} \\
            --talos_config ${talos_config} \\
            --run_date ${run_date} \\
            --output input.tsv
        """
}
