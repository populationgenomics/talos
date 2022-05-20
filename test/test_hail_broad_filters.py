"""
unit testing collection for the hail MT methods
"""

import os
from unittest.mock import MagicMock

import pytest
import hail as hl

from reanalysis.hail_filter_and_label import (
    filter_matrix_by_ac,
    filter_matrix_by_variant_attributes,
)


PWD = os.path.dirname(__file__)
INPUT = os.path.join(PWD, 'input')

hl_locus = hl.Locus(contig='chr1', position=1, reference_genome='GRCh38')


# I don't want a test fixture vcf with 70+ samples, so this is a pure Mock test
def test_filter_matrix_by_ac_small():
    """
    check the ac filter is not triggered
    """
    conf = {'min_samples_to_ac_filter': 10, 'ac_filter_percentage': 10}
    matrix_mock = MagicMock()
    # below threshold value
    matrix_mock.count_cols.return_value = 7
    mt = filter_matrix_by_ac(matrix_data=matrix_mock, config=conf)
    assert mt.filter_rows.call_count == 0


def test_filter_matrix_by_ac_large():
    """
    check the ac filter is triggered
    """
    conf = {'min_samples_to_ac_filter': 1, 'ac_threshold': 0.1}
    # above threshold value
    matrix_mock = MagicMock()
    matrix_mock.count_cols.return_value = 10
    matrix_mock.info.AC = 1
    matrix_mock.info.AN = 1
    matrix_mock.filter_rows.return_value = matrix_mock
    mt = filter_matrix_by_ac(matrix_data=matrix_mock, config=conf)
    assert mt.filter_rows.call_count == 1


@pytest.mark.parametrize(
    'filters,alleles,length',
    [
        (hl.empty_set(hl.tstr), hl.literal(['A', 'C']), 1),
        (hl.empty_set(hl.tstr), hl.literal(['A', '*']), 0),
        (hl.empty_set(hl.tstr), hl.literal(['A', 'C', 'G']), 0),
        (hl.literal({'fail'}), hl.literal(['A', 'C']), 0),
        (hl.literal({'VQSR'}), hl.literal(['A', 'C']), 0),
    ],
)
def test_filter_matrix_by_variant_attributes(filters, alleles, length, hail_matrix):
    """
    input 'values' are filters, & alleles

    :param filters:
    :param alleles:
    :param length:
    :param hail_matrix:
    :return:
    """
    # to add new alleles, we need to scrub alleles from the key fields
    hail_matrix = hail_matrix.key_rows_by('locus')
    anno_matrix = hail_matrix.annotate_rows(
        filters=filters,
        alleles=alleles,
    )

    anno_matrix = filter_matrix_by_variant_attributes(anno_matrix)
    assert anno_matrix.count_rows() == length
