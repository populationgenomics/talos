process MakeClinvarbitrationPm5 {
    container params.container

    publishDir params.processed_annotations, mode: 'copy'

    input:
        path annotated_snv

    def timestamp = new java.util.Date().format('yyyy-MM')

    output:
        path "clinvarbitration_${timestamp}.pm5.ht", emit: "ht"

    """
    python3 -m clinvarbitration.scripts.clinvar_by_codon \
        -i "${annotated_snv}" \
        -o "clinvarbitration_${timestamp}.pm5"
    rm clinvarbitration_${timestamp}.pm5.tsv
    """
}
