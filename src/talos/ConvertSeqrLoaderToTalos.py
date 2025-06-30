#!/usr/bin/python3

"""
This script intends to solve the problem of converting a Seqr Loader VEP annotated MT to a Talos-compatible format.
Original Issue: https://github.com/populationgenomics/talos/issues/516

the target schema to extract from the input MT are:

'transcript_consequences': array<struct {
    consequence: str,
    gene: str,
    transcript: str,
    biotype: str,
    amino_acid_change: str,
    dna_change: str,
    codon: int32,
    gene_id: str,
    am_class: str,
    am_pathogenicity: float64,
    mane_status: str,
    ensp: str,
    mane_id: str
},
'gene_ids': array<str>,

The VEP-annotated MT has a mostly equivalent set of fields, with some differences in naming and structure.

There are some fields with no explicit equivalent.

VEP-seqr-loader expected schema (taken from the issue):

transcript_consequences: array<struct {
    allele_num: int32,
    amino_acids: str,
    appris: str,
    biotype: str,
    canonical: int32,
    ccds: str,
    cdna_start: int32,
    cdna_end: int32,
    cds_end: int32,
    cds_start: int32,
    codons: str,
    consequence_terms: array<str>,
    distance: int32,
    domains: array<struct {
        db: str,
        name: str
    }>,
    exon: str,
    gene_id: str,
    gene_pheno: int32,
    gene_symbol: str,
    gene_symbol_source: str,
    hgnc_id: str,
    hgvsc: str,
    hgvsp: str,
    hgvs_offset: int32,
    impact: str,
    intron: str,
    lof: str,
    lof_flags: str,
    lof_filter: str,
    lof_info: str,
    minimised: int32,
    polyphen_prediction: str,
    polyphen_score: float64,
    protein_end: int32,
    protein_start: int32,
    protein_id: str,
    sift_prediction: str,
    sift_score: float64,
    strand: int32,
    swissprot: str,
    transcript_id: str,
    trembl: str,
    tsl: int32,
    uniparc: str,
    variant_allele: str
}>,
"""

from argparse import ArgumentParser

import hail as hl


# lookup for converting VEP <-> BCFtools consequence terms
VEP_TO_CSQ = {
    'frameshift_variant': 'frameshift',
    'splice_acceptor_variant': 'splice_acceptor',
    'splice_donor_variant': 'splice_donor',
}


def main(input_mt: str, output_mt: str):
    """Reformat annotations on a Seqr Loader VEP annotated MT to Talos-compatible format."""

    # replace this initialization with your Hail context initialization of choice - batch backend? local?
    hl.init()

    # Load the input MT
    mt = hl.read_matrix_table(input_mt)

    consequence_translation = hl.literal(VEP_TO_CSQ)

    # test, to delete
    mt = mt.annotate_rows(
        transcript_consequences=mt.transcript_consequences.map(
            lambda x: hl.struct(
                consequence=hl.delimit(
                    hl.map(
                        lambda term: consequence_translation.get(term, term),  # translate VEP terms to CSQ terms
                        x.consequence.split('&'),
                    ),
                    '&',
                )
            ),
        ),
        gene_ids=mt.gene_ids,  # this is already an array of strings
    )

    # region: VEP annotations
    mt = mt.annotate_rows(
        transcript_consequences=mt.vep.transcript_consequences.map(
            lambda x: hl.struct(
                consequence=hl.delimit(
                    hl.map(
                        lambda term: consequence_translation.get(term, term),  # translate VEP terms to CSQ terms
                        x.consequence_terms,
                    ),
                    '&',
                ),
                gene=x.gene_id,
                transcript=x.transcript_id,
                biotype=x.biotype,
                amino_acid_change=x.hgvsp,  # HGVSP is the annotation we want, BCFTools doesn't supply this
                dna_change=x.hgvsc,  # HGVSC is preferable here, BCFTools is again inferior
                codon=x.protein_start,  # a better replacement would be x.codons.split('/')[0] if codons are present
                # these are both VEP-annotations in VEP 110
                # https://github.com/populationgenomics/production-pipelines/blob/main/cpg_workflows/query_modules/vep.py#L188-L189
                am_class=None,
                am_pathogenicity=None,
                # https://github.com/populationgenomics/production-pipelines/blob/main/cpg_workflows/query_modules/vep.py#L116-L117
                # For us these are status (Select/Plus Clinical) and ID (NMID)
                mane_status=None,
                mane_id=None,
                ensp=x.protein_id,
            ),
        ),
        gene_ids=mt.vep.gene_ids,  # this is already an array of strings
    )
    # endregion

    # region: gnomAD
    # pull out gnomAD annotations? expected struct here is:
    _exp_struct_schema = """
    'gnomad': struct {
        gnomad_AC: int32,
        gnomad_AF: float64,
        gnomad_AC_XY: int32,
        gnomad_HomAlt: int32
    }
    """
    # we used to get this from a seqr_loader companion table, now we get it using echtvar. Unsure of where this sits
    # in your MTs, only the transcript consequences were specified in the issue.
    # endregion

    # drop the VEP annotations
    mt.drop('vep')

    # region: extra actions prior to export
    # additional considerations before saving a new object:
    # 1. removing unwanted samples: https://github.com/populationgenomics/talos/blob/1f7f50a31733188158c1fe3bf90236c664cbd8e7/src/talos/RunHailFiltering.py#L865-L908
    # 2. removing unwanted rows: https://github.com/populationgenomics/talos/blob/main/src/talos/cpg_internal_scripts/extract_fragmented_vcf_from_mt.py#L72-L77
    # endregion

    # Save the modified MT to the output path
    mt.write(output_mt, overwrite=True)


if __name__ == '__main__':
    parser = ArgumentParser(description='Convert a Seqr Loader VEP annotated MT to Talos-compatible format.')
    parser.add_argument('input_mt', type=str, help='Path to the input MT file')
    parser.add_argument('output_mt', type=str, help='Path to the output MT file')
    args = parser.parse_args()

    main(input_mt=args.input_mt, output_mt=args.output_mt)
