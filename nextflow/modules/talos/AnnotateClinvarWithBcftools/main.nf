process AnnotateClinvarWithBcftools {
    container params.container

    publishDir params.processed_annotations

    input:
        path vcf
        path ref_fa
        path gff3

    output:
        path "clinvarbitration.annotated.tsv", emit: tsv

    """
    bcftools csq \
        --force \
        -f "${ref_fa}" \
        -g "${gff3}" \
        --unify-chr-names 'chr,-,chr' \
        "${vcf}" |
    bcftools +split-vep \
        -d \
        -s :missense \
        -f "%transcript\t%amino_acid_change\t%allele_id\t%gold_stars\n" \
        - \
        > clinvarbitration.annotated.tsv
    """
}
