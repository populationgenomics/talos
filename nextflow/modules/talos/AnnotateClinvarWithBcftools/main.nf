process AnnotateClinvarWithBcftools {
    container params.container

    publishDir params.output_dir

    input:
        path vcf
        path ref_fa
        path gff3

    output:
        path "clinvarbitration.annotated.tsv", emit: tsv

    """
    bcftools csq \
        -f "${ref_fa}" \
        -g "${gff3}" \
        "${vcf}" |
    bcftools +split-vep \
        -d \
        -s :missense \
        -f "%transcript\t%amino_acid_change\t%allele_id\t%gold_stars\n" \
        - \
        > clinvarbitration.annotated.tsv
    """
}
