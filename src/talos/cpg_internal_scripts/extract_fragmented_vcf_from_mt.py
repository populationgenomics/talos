#!/usr/bin/env python3

"""
Takes a MatrixTable and two output paths
Writes two representations - region filtered MatrixTable, and the same as a sites-only VCF representation
VCF is exported per-partition, with a separate header file, to be concatenated later

This script also takes a BED file as input; output contains only the variants that overlap with the BED file

All existing INFO fields are dropped, and replaced with just the callset AC / AN / AF
"""

import argparse

import loguru

import hail as hl

from cpg_utils import config, hail_batch

# Update the VQSR header elements so 3rd party tools can read the VCF reliably
VQSR_FILTERS = {
    'filter': {
        'VQSRTrancheINDEL99.00to99.50': {
            'Description': 'Truth sensitivity tranche level for INDEL model at VQS Lod: -1.4652 <= x < -0.6489',
        },
        'VQSRTrancheINDEL99.50to99.90': {
            'Description': 'Truth sensitivity tranche level for INDEL model at VQS Lod: -8.3914 <= x < -1.4652',
        },
        'VQSRTrancheINDEL99.90to99.95': {
            'Description': 'Truth sensitivity tranche level for INDEL model at VQS Lod: -20.9224 <= x < -8.3914',
        },
        'VQSRTrancheINDEL99.95to100.00': {
            'Description': 'Truth sensitivity tranche level for INDEL model at VQS Lod: -39995.8675 <= x < -20.9224',
        },
        'VQSRTrancheINDEL99.95to100.00+': {
            'Description': 'Truth sensitivity tranche level for INDEL model at VQS Lod < -39995.8675',
        },
        'VQSRTrancheSNP99.00to99.90+': {
            'Description': 'Truth sensitivity tranche level for SNP model at VQS Lod < -10.0',
        },
        'VQSRTrancheSNP99.90to100.00': {
            'Description': 'Truth sensitivity tranche level for SNP model at VQS Lod: -10.0 <= x < -4.37',
        },
        'VQSRTrancheSNP99.90to100.00+': {
            'Description': 'Truth sensitivity tranche level for SNP model at VQS Lod < -10.0',
        },
        'MONOALLELIC': {
            'Description': 'Variant is monoallelic?',
        },
    },
}


def main(
    mt_path: str,
    output_mt: str,
    output_sites_only: str,
    bed: str | None,
) -> None:
    """

    Args:
        mt_path (str):
        output_mt (str): write region-filtered MatrixTable, stripped of INFO fields
        output_sites_only (str): write a per-partition sites-only VCF directory to this location
        bed (str): Region BED file
    """

    hail_batch.init_batch()

    # read the dense MT and obtain the sites-only HT
    mt = hl.read_matrix_table(mt_path)

    if bed:
        # remote-read of the BED file, skipping any contigs not in the reference genome
        # the Ensembl data wasn't filtered to remove non-standard contigs
        limited_region = hl.import_bed(bed, reference_genome=config.genome_build(), skip_invalid_intervals=True)
        # filter to overlaps with the BED file
        mt = mt.filter_rows(hl.is_defined(limited_region[mt.locus]))

    # replace the existing INFO block to just have AC/AN/AF - no other carry-over. Allow for this to be missing.
    if 'AF' not in mt.info:
        mt = hl.variant_qc(mt)
        mt = mt.annotate_rows(
            info=mt.info.annotate(
                AF=[mt.variant_qc.AF[1]],
                AN=mt.variant_qc.AN,
                AC=[mt.variant_qc.AC[1]],
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

    mt.describe()

    mt.write(output_mt, overwrite=True)

    # now read that location for speed, and write the sites-only VCF
    # keep partitions consistent
    sites_only_ht = hl.read_matrix_table(output_mt).rows()

    loguru.logger.info('Writing sites-only VCF in fragments, header-per-shard')
    hl.export_vcf(
        sites_only_ht,
        output_sites_only,
        tabix=True,
        parallel='separate_header',
        metadata=VQSR_FILTERS,
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input',
        help='Path to the input MT',
        required=True,
    )
    parser.add_argument(
        '--output_mt',
        help='Path to write the resulting MatrixTable',
        required=True,
    )
    parser.add_argument(
        '--output_sites_only',
        help='Specify an output path for a sites-only VCF, or None',
    )
    parser.add_argument(
        '--bed',
        help='Region BED file',
    )
    args = parser.parse_args()
    main(mt_path=args.input, output_mt=args.output_mt, output_sites_only=args.output_sites_only, bed=args.bed)
