#!/usr/bin/env python3

"""
Takes the full un-annotated joint-callset VCF,
A HailTable of the formatted annotations,
Integrates the two, writing a MatrixTable representation of the fully annotated VCF
"""

import logging
from argparse import ArgumentParser
from pathlib import Path

import hail as hl


def cli_main():
    """
    take an input VCF and an output MT path
    also supply the alpha_missense table created by parse_amissense_into_ht.py
    """

    parser = ArgumentParser(description='Takes a HT of annotations, and a callset VCF, and combines into a MT')
    parser.add_argument('--input', help='Path to the annotated sites-only VCF', required=True)
    parser.add_argument('--annotations', help='Path to the annotated sites-only VCF', required=True)
    parser.add_argument('--output', help='output Table path, must have a ".ht" extension', required=True)
    args = parser.parse_args()

    assert args.output.endswith('.mt'), 'Output path must end in .mt'

    # check specifically for a SUCCESS file, marking a completed hail write
    # will fail if we accidentally pass the compressed Tarball path
    assert (Path(args.annotations) / '_SUCCESS').exists(), 'Annotations Table does not exist'

    main(
        vcf_path=args.input,
        output_path=args.output,
        annotations=args.annotations,
    )


def main(
    vcf_path: str,
    output_path: str,
    annotations: str,
):
    """
    Takes a Hail-Table of annotations, a joint-called VCF, reads the VCF as a MatrixTable and hops the annotations over

    I'm initially trying this without a checkpoint, though I expect that not to scale well

    Args:
        vcf_path (str): path to the full callset VCF
        output_path (str): path to write the resulting MatrixTable to
        annotations (str): path to a Hail Table containing annotations
    """

    hl.default_reference('GRCh38')

    # read the VCF into a MatrixTable
    mt = hl.import_vcf(
        vcf_path,
        array_elements_required=False,
        force_bgz=True,
    )

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

    mt.write(output_path)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    cli_main()
