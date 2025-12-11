"""
This script takes an initial VCF, and converts it to a MT
The previous process had some slow VCF operations, then converted the VCF to a MT anyway
"""

from argparse import ArgumentParser

from loguru import logger

import hail as hl


def main(input_path: str, bed: str, output_path: str) -> None:
    hl.context.init_spark(master='local[*]', default_reference='GRCh38')
    mt = hl.import_vcf(
        input_path,
        array_elements_required=False,
        force_bgz=True,
    ).checkpoint(f'{output_path}._temp', _read_if_exists=True)

    # allow for IGG-related shennanigans
    mt = mt.annotate_entries(
        AD=hl.if_else(
            # if this is a high quality HomRef
            (mt.GT.is_hom_ref()) & (hl.len(mt.AD) == 1),
            hl.if_else(
                # If HomRef is populated
                hl.is_defined(mt.AD),
                mt.AD.append(0),
                mt.AD,
            ),
            mt.AD,
        ),
    )

    logger.info(f'Reading the bed file {bed} to use in region-filtering {input_path}')
    # read the BED file, skipping any contigs not in the reference genome
    bed_region = hl.import_bed(bed, skip_invalid_intervals=True)

    # filter to overlaps with the BED file
    mt = mt.filter_rows(hl.is_defined(bed_region[mt.locus]))

    logger.info('Splitting multiallelic variants in input data')
    mt = hl.split_multi_hts(mt)

    # replace the existing INFO block to just have AC/AN/AF - no other carry-over. Allow for this to be missing.
    if 'AF' not in mt.info:
        mt = hl.variant_qc(mt)
        mt = mt.annotate_rows(
            info=hl.struct(
                AF=mt.variant_qc.AF,
                AN=mt.variant_qc.AN,
                AC=mt.variant_qc.AC,
            ),
            filters=hl.empty_set(hl.tstr),
        )
        mt = mt.drop('variant_qc')

    mt = mt.select_rows(
        info=hl.struct(
            AF=mt.info.AF,
            AN=mt.info.AN,
            AC=mt.info.AC,
        ),
        rsid=mt.rsid,
        filters=mt.filters,
    )

    mt.write(output_path)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', help='VCF file to convert into a MT')
    parser.add_argument('--regions', help='Bed file to filter by')
    parser.add_argument('--output', help='Output MT')
    args = parser.parse_args()
    main(input_path=args.input, bed=args.regions, output_path=args.output)
