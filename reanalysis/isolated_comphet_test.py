"""
single script to check the compound-het failure on minimal data
"""

from typing import Dict, List
from itertools import permutations
import json
import logging
import sys
from argparse import ArgumentParser

import hail as hl

from cloudpathlib import AnyPath
from cpg_utils.hail import init_batch


def transform_variant_string(locus_details: hl.Struct) -> str:
    """
    takes an object
    Struct(
        locus=Locus(
            contig='chr1',
            position=10,
            reference_genome='GRCh38'
        ),
        alleles=['GC', 'G'],
        category_4_only=0
    )

    transform into simplified 1-10-GC-G
    drop the category_4_only attribute
    :param locus_details:
    :return:
    """
    return '-'.join(
        [
            locus_details.locus.contig.replace('chr', ''),
            str(locus_details.locus.position),
            *locus_details.alleles,
        ]
    )


def extract_comp_het_details(
    matrix: hl.MatrixTable,
) -> Dict[str, Dict[str, Dict[str, List[str]]]]:
    """
    takes the matrix table, and finds compound-hets per sample
    based on the gene name only

    return format is a nested dictionary:
    Sample:
        Gene:
            Var1: [Var2, VarN],
            ..
        ..
    ..

    :param matrix:
    """

    logging.info('Extracting out the compound-het variant pairs')

    # set a new group of values as the key, so that we can collect on them easily
    ch_matrix = matrix.key_rows_by(matrix.locus, matrix.alleles, matrix.category_4_only)
    ch_matrix = ch_matrix.annotate_cols(
        hets=hl.agg.group_by(
            ch_matrix.info.gene_id,
            hl.agg.filter(ch_matrix.GT.is_het(), hl.agg.collect(ch_matrix.row_key)),
        )
    )

    # extract those possible compound het pairs out as a non-Hail structure
    compound_hets = {}

    logging.info('Collecting all variant pairs')

    # iterate over the hail table rows
    # find all variant pair permutations which aren't both class 4
    for row in ch_matrix.select_cols('hets').col.collect():

        # prepare a summary dict for this sample
        sample_dict = {}

        # iterate over all the `gene: [var1, var2]` structures
        for gene, variants in dict(row.hets).items():

            # assess each possible variant pairing
            for var1, var2 in permutations(variants, 2):

                # skip if both are class 4 only - not valuable pairing
                if var1.category_4_only == 1 and var2.category_4_only == 1:
                    continue

                # pair the string transformation
                sample_dict.setdefault(gene, {}).setdefault(
                    transform_variant_string(var1), []
                ).append(transform_variant_string(var2))

        # if we found comp hets, add the content for this sample
        if len(sample_dict) > 0:
            compound_hets[row.s] = sample_dict

    return compound_hets


def main(mt_in: str, out_json: str):
    """

    :param mt_in:
    :param out_json:
    :return:
    """

    # initiate Hail with upgraded driver spec.
    init_batch(
        driver_cores=8,
        driver_memory='highmem',
    )

    # load in the MT
    matrix = hl.read_matrix_table(mt_in)

    # parse out the compound het details (after pulling gene_id above)
    comp_het_details = extract_comp_het_details(matrix=matrix)

    # transform the vcf output path into a json path
    out_json = f'{out_json.split(".", maxsplit=1)[0]}.json'

    # and write the comp-het JSON file
    serialised_obj = json.dumps(comp_het_details, indent=True, default=str)
    AnyPath(out_json).write_text(serialised_obj)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )

    parser = ArgumentParser()
    parser.add_argument(
        '--mt_input',
        required=True,
        help='path to the matrix table to ingest',
    )
    parser.add_argument(
        '--out_json', type=str, required=True, help='VCF path to export results'
    )
    args = parser.parse_args()
    main(
        mt_in=args.mt_input,
        out_json=args.out_json,
    )
