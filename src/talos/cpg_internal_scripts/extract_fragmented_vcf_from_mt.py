#!/usr/bin/env python3

"""
Takes a MatrixTable, file containing SG IDs, and an output path

Reads the MT, filters to SG IDs in the file

Writes the VCF out in shards, each fragment containing a full header-per-shard

All existing INFO fields are dropped, and replaced with just the callset AC / AN / AF
"""

import argparse

import loguru

import hail as hl

from cpg_utils import config, hail_batch, to_path

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
        'LowQual': {
            'Description': 'Low Quality Variant...',
        },
    },
}


def filter_mt_to_sgids(
    mt: hl.MatrixTable,
    sgid_file: str,
) -> hl.MatrixTable:
    """Read the full MatrixTable, and subset to a collection of SG IDs in a text file."""

    # read SG IDs from a file
    id_file = to_path(sgid_file)
    if not id_file.exists():
        raise ValueError(f'Sample ID file {id_file} does not exist')

    id_list = {each.strip() for each in id_file.read_text().splitlines()}

    mt_sample_ids = set(mt.s.collect())

    # fail if there
    if sample_ids_not_in_mt := id_list - mt_sample_ids:
        raise ValueError(
            f'Found {len(sample_ids_not_in_mt)}/{len(id_list)} IDs in the requested subset not in the callset.\n'
            f"IDs that aren't in the callset: {sample_ids_not_in_mt}\n"
            f'All callset sample IDs: {mt_sample_ids}',
        )

    loguru.logger.info(f'Found {len(mt_sample_ids)} samples in mt, subsetting to {len(id_list)} samples.')

    mt = mt.filter_cols(hl.literal(id_list).contains(mt.s))
    return mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))


def main(
    mt_path: str,
    sg_id_file: str,
    output: str,
) -> None:
    """

    Args:
        mt_path (str):
        sg_id_file (str): file containing SG IDs
        output (str): write VCFs, stripped of INFO fields and annotations
    """
    hail_batch.init_batch(
        worker_memory=config.config_retrieve(['combiner', 'worker_memory'], 'highmem'),
        worker_cores=config.config_retrieve(['combiner', 'worker_cores'], 2),
        driver_memory=config.config_retrieve(['combiner', 'driver_memory'], 'highmem'),
        driver_cores=config.config_retrieve(['combiner', 'driver_cores'], 2),
    )

    # read the dense MT and obtain the sites-only HT
    mt = hl.read_matrix_table(mt_path)

    # filter the MT to a specific set of SG IDs
    mt = filter_mt_to_sgids(mt, sg_id_file)

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

    # drop all annotations
    mt = mt.select_rows(
        info=hl.struct(
            AF=mt.info.AF,
            AN=mt.info.AN,
            AC=mt.info.AC,
        ),
        rsid=mt.rsid,
        filters=mt.filters,
    )

    # determine how many fragments to generate - the thousands we have by default will be pretty excessive
    if fragments := config.config_retrieve(['workflow', 'vcf_fragments'], False):
        mt = mt.repartition(fragments)

    loguru.logger.info('Writing sites-only VCF in fragments, header-per-shard')
    hl.export_vcf(
        mt,
        output,
        tabix=True,
        parallel='header_per_shard',
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
        '--sgs',
        help='Path to the SG ID file',
        required=True,
    )
    parser.add_argument(
        '--output',
        help='Path to write the resulting VCF',
        required=True,
    )
    args = parser.parse_args()
    main(mt_path=args.input, sg_id_file=args.sgs, output=args.output)
