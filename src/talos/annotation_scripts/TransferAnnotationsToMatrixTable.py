#!/usr/bin/env python3

"""
Takes the full un-annotated joint-callset VCF,
A HailTable of the formatted annotations,
Integrates the two, writing a MatrixTable representation of the fully annotated VCF
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

    parser = ArgumentParser(description='Takes a HT of annotations, and a callset VCF, and combines into a MT')
    parser.add_argument('--input', help='Path to the MatrixTable', required=True)
    parser.add_argument('--annotations', help='Path to the annotated sites-only VCF', required=True)
    parser.add_argument('--output', help='output Table path, must have a ".ht" extension', required=True)
    parser.add_argument('--backend', help='type of backend to use', default='local')
    args = parser.parse_args()

    main(
        input_path=args.input,
        output_path=args.output,
        annotations=args.annotations,
        backend=args.backend,
    )


def main(
    input_path: str,
    output_path: str,
    annotations: str,
    backend: str,
):
    """
    Takes a Hail-Table of annotations, a joint-called VCF, reads the VCF as a MatrixTable and hops the annotations over

    I'm initially trying this without a checkpoint, though I expect that not to scale well

    Args:
        input_path (str): path to the full callset VCF, or MT of the same
        output_path (str): path to write the resulting MatrixTable to
        annotations (str): path to a Hail Table containing annotations
        backend (str): which backend type to use
    """

    if backend == 'local':
        logger.info('Using local backend for Hail')
        hl.context.init_spark(master='local[*]', default_reference='GRCh38', quiet=True)
    else:
        logger.info('Using Batch backend for Hail')
        init_batch()

    # read the VCF into a MatrixTable, or read the existing MatrixTable
    if input_path.endswith('.mt'):
        logger.info(f'Reading MatrixTable from {input_path!r}')
        mt = hl.read_matrix_table(input_path)
    elif input_path.endswith(('.vcf', '.vcf.gz', '.vcf.bgz')):
        logger.info(f'Reading VCF from {input_path!r}')
        # import the VCF, and set the array_elements_required to False
        # this is because the VCF may not have all fields present in all rows
        # which is fine for our purposes
        mt = hl.import_vcf(
            input_path,
            array_elements_required=False,
            force_bgz=True,
        )
    else:
        raise ValueError(f'Input path must be a VCF or MatrixTable, got {input_path!r}')

    # read the annotations into a Table
    ht = hl.read_table(annotations)

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

    mt.write(output_path, overwrite=True)


if __name__ == '__main__':
    cli_main()
