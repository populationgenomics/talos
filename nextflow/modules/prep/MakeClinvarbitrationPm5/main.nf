process MakeClinvarbitrationPm5 {
    container params.container

    input:
        path annotated_snv
        val timestamp

    output:
        path "clinvarbitration_${timestamp}.pm5.ht"

    script:
    """
    python3 -m clinvarbitration.scripts.clinvar_by_codon \
        -i "${annotated_snv}" \
        -o "clinvarbitration_${timestamp}.pm5"
    rm clinvarbitration_${timestamp}.pm5.tsv
    """
}
