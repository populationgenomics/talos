
process ConvertSpliceVarDb {
    container params.container

    publishDir params.cohort_output_dir, mode: 'copy'

    input:
        path svdb

    output:
        path "${params.cohort}_svdb.ht.tar.gz"

    """
    python -m talos.convert_splice_var_db \
        --input ${svdb} \
        --output ${params.cohort}_svdb.ht
    tar --no-xattrs -czf ${params.cohort}_svdb.ht.tar.gz ${params.cohort}_svdb.ht
    """
}
