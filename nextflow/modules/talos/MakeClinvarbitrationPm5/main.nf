process MakeClinvarbitrationPm5 {
    container params.container

    publishDir params.processed_annotations, mode: 'copy'

    input:
        path annotated_snv

    output:
        path "clinvarbitration.pm5.ht", emit: "ht"

    """
    python3 -m clinvarbitration.scripts.clinvar_by_codon \
        -i "${annotated_snv}" \
        -o "clinvarbitration.pm5"
    rm clinvarbitration.pm5.tsv
    """
}
