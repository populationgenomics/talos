#!/usr/bin/env python3

"""
Takes the MatrixTable representing the full callset, and one or more Hail Tables containing annotations
Integrates the two, writing a MatrixTable representation of the fully annotated joint call
"""

from argparse import ArgumentParser

from loguru import logger

import hail as hl

from cpg_utils.hail_batch import init_batch


def cli_main():
    """
    take an input VCF and an output MT path
    also supply the alpha_missense table created by parse_amissense_into_ht.py
    """

    parser = ArgumentParser(description='Takes a HT of annotations, and a callset MT, and combines both')
    parser.add_argument('--input', help='Path to the input MatrixTable', required=True)
    parser.add_argument(
        '--annotations',
        help='Path to one or more annotation tables',
        required=True,
        nargs='+',
    )
    parser.add_argument('--output', help='output Table path, must have a ".ht" extension', required=True)
    parser.add_argument('--backend', help='type of backend to use', default='local')
    parser.add_argument('--regions', help='Optional, BED file to write reduced output', default=None)
    args = parser.parse_args()

    main(
        input_path=args.input,
        output_path=args.output,
        annotations=args.annotations,
        backend=args.backend,
        regions=args.regions,
    )


def main(
    input_path: str,
    output_path: str,
    annotations: list[str],
    backend: str,
    regions: str | None = None,
):
    """
    Takes a Hail-Table of annotations, a joint-called VCF, reads the VCF as a MatrixTable and hops the annotations over

    I'm initially trying this without a checkpoint, though I expect that not to scale well

    Args:
        input_path (str): path to the full callset VCF, or MT of the same
        output_path (str): path to write the resulting MatrixTable to
        annotations (list[str]): path to one or more Hail Tables containing annotations
        backend (str): which backend type to use
        regions (str | None): optional, a BED file to limit the written data to only genic regions (space saving)
    """

    if not annotations:
        raise ValueError('No annotations HTs provided')

    if backend == 'local':
        logger.info('Using local backend for Hail')
        hl.context.init_spark(master='local[*]', default_reference='GRCh38', quiet=True)
    else:
        logger.info('Using Batch backend for Hail')
        init_batch()

    # read the VCF into a MatrixTable, or read the existing MatrixTable
    logger.info(f'Reading MatrixTable from {input_path!r}')
    mt = hl.read_matrix_table(input_path)

    # read each annotation table individually, then obtain the union of all tables
    annotation_tables = [hl.read_table(each_table) for each_table in annotations]
    if len(annotation_tables) == 1:
        ht = annotation_tables[0]
    else:
        ht = annotation_tables[0].union(*annotation_tables[1:])

    # syntax sweeter for later on
    matched_annotations = ht[mt.row_key]

    # a couple of lines commented of to make this as easy as possible to adopt
    # talos doesn't make use of these annotations yet
    # ruff: noqa: ERA001
    mt = mt.annotate_rows(
        gnomad=hl.struct(
            gnomad_AC=matched_annotations.info.gnomad_AC_joint,
            gnomad_AF=matched_annotations.info.gnomad_AF_joint,
            # gnomad_AN=matched_annotations.info.gnomad_AN_joint,
            gnomad_AC_XY=matched_annotations.info.gnomad_AC_joint_XY,
            # gnomad_AF_XY=matched_annotations.info.gnomad_AF_joint_XY,
            # gnomad_FAF=matched_annotations.info.gnomad_faf_95_joint,
            gnomad_HomAlt=matched_annotations.info.gnomad_HomAlt_joint,
        ),
        transcript_consequences=matched_annotations.transcript_consequences,
        gene_ids=matched_annotations.gene_ids,
    )

    mt.describe()

    if regions:
        logger.info(f'Reading the bed file {regions} to use in region-filtering {input_path}')
        # read of the BED file, skipping any contigs not in the reference genome
        bed_region = hl.import_bed(regions, skip_invalid_intervals=True)

        # filter to overlaps with the BED file
        mt = mt.filter_rows(hl.is_defined(bed_region[mt.locus]))

    mt.write(output_path, overwrite=True)


if __name__ == '__main__':
    cli_main()
