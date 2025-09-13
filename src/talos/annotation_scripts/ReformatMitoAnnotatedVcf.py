#!/usr/bin/env python3

"""
This is an adapter process to take a sites-only VCF annotated with gnomAD frequencies and BCFtools csq consequences, and
re-arrange it into a HailTable for use with the Talos pipeline.

This process combines the AF/CSQs already applied with the MANE transcript/protein names, and AlphaMissense annotations
"""

import json
from argparse import ArgumentParser
from collections import defaultdict

from cpg_utils.hail_batch import init_batch
from loguru import logger

import hail as hl

from talos.annotation_scripts.ReformatAnnotatedVcfIntoHailTable import (
    extract_and_split_csq_string,
    csq_strings_into_hail_structs,
    MISSING_FLOAT,
    MISSING_INT,
    MISSING_STRING,
)


def csq_strings_into_hail_structs(csq_strings: list[str], mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Take the list of BCSQ strings, split the CSQ annotation and re-organise as a hl struct

    Args:
        csq_strings (list[str]): a list of strings, each representing a CSQ entry
        mt (hl.Table): the Table to annotate

    Returns:
        a Table with the BCSQ annotations re-arranged
    """

    # get the BCSQ contents as a list of lists of strings, per variant
    split_csqs = mt.info.BCSQ.map(lambda csq_entry: csq_entry.split('\|'))  # noqa: W605

    # this looks pretty hideous, bear with me
    # if BCFtools csq doesn't have a consequence annotation, it will truncate the pipe-delimited string
    # this is fine sometimes, but not when we're building a schema here
    # when we find truncated BCSQ strings, we need to add dummy values to the end of the array
    split_csqs = split_csqs.map(
        lambda x: hl.if_else(
            # if there were only 4 values, add 3 missing Strings
            hl.len(x) == 4,
            x.extend([MISSING_STRING, MISSING_STRING, MISSING_STRING]),
            hl.if_else(
                # 5 values... add 2 missing Strings
                hl.len(x) == 5,
                x.extend([MISSING_STRING, MISSING_STRING]),
                hl.if_else(
                    hl.len(x) == 6,
                    x.extend([MISSING_STRING]),
                    x,
                ),
            ),
        ),
    )

    # transform the CSQ string arrays into structs using the header names
    # Consequence | gene | transcript | biotype | strand | amino_acid_change | dna_change
    mt = mt.annotate_rows(
        transcript_consequences=split_csqs.map(
            lambda x: hl.struct(
                **{csq_strings[n]: x[n] for n in range(len(csq_strings)) if csq_strings[n] != 'strand'},
            ),
        ),
    )

    return mt.annotate_rows(
        # amino_acid_change can be absent, or in the form of "123P" or "123P-124F"
        # we use this number when matching to the codons of missense variants, to find codon of the reference pos.
        transcript_consequences=hl.map(
            lambda x: x.annotate(
                codon=hl.if_else(
                    x.amino_acid_change == MISSING_STRING,
                    hl.missing(hl.tint32),
                    hl.if_else(
                        x.amino_acid_change.matches('^([0-9]+).*$'),
                        hl.int32(x.amino_acid_change.replace('^([0-9]+).+', '$1')),
                        hl.missing(hl.tint32),
                    ),
                ),
                am_class=MISSING_STRING,
                am_pathogenicity=MISSING_FLOAT,
                mane_status=MISSING_STRING,
                ensp=MISSING_STRING,
                mane_id=MISSING_STRING,
            ),
            mt.transcript_consequences,
        ),
    )


def cli_main():
    """
    take an input VCF and an output MT path
    also supply the alpha_missense table created by parse_amissense_into_ht.py
    """

    parser = ArgumentParser(description='Takes a BCSQ annotated VCF and makes it a mt')
    parser.add_argument('--input', help='Path to the annotated sites-only VCF', required=True)
    parser.add_argument('--output', help='output Table path, must have a ".mt" extension', required=True)
    parser.add_argument('--batch', help='flag to use the batch hail backend', action='store_true')
    args = parser.parse_args()

    main(
        vcf_path=args.input,
        output_path=args.output,
        batch=args.batch,
    )


def main(
    vcf_path: str,
    output_path: str,
    batch: bool,
):
    """
    Takes a VEP-annotated Mito VCF, reorganises into a Talos-compatible MatrixTable

    Args:
        vcf_path (str): path to the annotated sites-only VCF
        output_path (str): path to write the resulting Hail Table to, must
        batch (bool): whether or not to batch the Hail Batch backend
    """

    if batch:
        logger.info('Using Hail Batch backend')
        init_batch()
    else:
        logger.info('Using Hail Local backend')
        hl.context.init_spark(master='local[*]', default_reference='GRCh38', quiet=True)

    # pull and split the CSQ header line
    csq_fields = extract_and_split_csq_string(vcf_path=vcf_path)

    # read the VCF into a MatrixTable
    mt = hl.import_vcf(vcf_path, array_elements_required=False, force_bgz=True)

    # re-shuffle the BCSQ elements
    mt = csq_strings_into_hail_structs(csq_fields, mt)

    # todo no ENSG IDs currently available - investigate MANE MT genes
    mt = mt.annotate_rows(
        transcript_consequences=hl.map(
            lambda x: x.annotate(
                gene_id=x.gene,
            ),
            mt.transcript_consequences,
        ),
        gnomad=hl.struct(
            gnomad_AC=MISSING_INT,
            gnomad_AF=MISSING_FLOAT,
            gnomad_AC_XY=MISSING_INT,
            gnomad_HomAlt=MISSING_INT,
        ),
    )

    # get a hold of the geneIds - use some aggregation
    mt = mt.annotate_rows(gene_ids=hl.set(mt.transcript_consequences.map(lambda c: c.gene_id)))

    # drop the BCSQ field
    mt = mt.annotate_rows(info=mt.info.drop('BCSQ'))

    mt.describe()

    mt.write(output_path, overwrite=True)


if __name__ == '__main__':
    cli_main()
