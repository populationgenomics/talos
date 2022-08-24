#!/usr/bin/env python3

"""
Prepare test data for annotation tests: subset the annotation refernece hail tables to the ROI
"""

import hail as hl
from cpg_utils.hail_batch import init_batch

init_batch()

chrom = 'chr20'
pos = '5111495-5111607'

ht = hl.read_table('gs://cpg-reference/seqr/v0-1/clinvar.GRCh38.ht')
ht = hl.filter_intervals(ht, [hl.parse_locus_interval(f'{chrom}:{pos}')])
ht.write_table(f'gs://cpg-fewgenomes-test/unittest/input/clinvar.GRCh38-{chrom}-{pos}.ht')
print('Done writing clinvar')

ht = hl.read_table(
    'gs://cpg-reference/seqr/v0-1/combined_reference_data_grch38-2.0.4.ht'
)
ht = hl.filter_intervals(ht, [hl.parse_locus_interval(f'{chrom}:{pos}')])
ht.write_table(
    f'gs://cpg-fewgenomes-test/unittest/input/combined_reference_data_grch38-2.0.4-{chrom}-{pos}.ht'
)
print('Done writing combined_reference_data_grch38')
