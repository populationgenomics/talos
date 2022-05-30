"""
unit testing collection for the hail MT methods
"""


from unittest.mock import MagicMock

import pytest
import hail as hl

from reanalysis.hail_filter_and_label import (
    filter_matrix_by_ac,
    filter_on_quality_flags,
    filter_to_well_normalised,
)


hl_locus = hl.Locus(contig='chr1', position=1, reference_genome='GRCh38')


def test_filter_matrix_by_ac_large():
    """
    check the ac filter is triggered
    """
    # above threshold value
    matrix_mock = MagicMock()
    matrix_mock.count_cols.return_value = 10
    matrix_mock.info.AC = 1
    matrix_mock.info.AN = 1
    matrix_mock.filter_rows.return_value = matrix_mock
    mt = filter_matrix_by_ac(matrix=matrix_mock)
    assert mt.filter_rows.call_count == 1


@pytest.mark.parametrize(
    'filters,length',
    [
        (hl.empty_set(hl.tstr), 1),
        (hl.literal({'fail'}), 0),
        (hl.literal({'VQSR'}), 0),
    ],
)
def test_filter_on_quality_flags(filters, length, hail_matrix):
    """
    annotate filters and run tests

    :param filters:
    :param length:
    :param hail_matrix:
    :return:
    """
    # to add new alleles, we need to scrub alleles from the key fields
    hail_matrix = hail_matrix.key_rows_by('locus')
    anno_matrix = hail_matrix.annotate_rows(filters=filters)

    assert filter_on_quality_flags(anno_matrix).count_rows() == length


@pytest.mark.parametrize(
    'alleles,length',
    [
        (hl.literal(['A', 'C']), 1),
        (hl.literal(['A', '*']), 0),
        (hl.literal(['A', 'C', 'G']), 0),
    ],
)
def test_filter_to_well_normalised(alleles, length, hail_matrix):
    """
    checks the allele-level tests

    :param alleles:
    :param length:
    :param hail_matrix:
    :return:
    """
    # to add new alleles, we need to scrub alleles from the key fields
    hail_matrix = hail_matrix.key_rows_by('locus')
    anno_matrix = hail_matrix.annotate_rows(alleles=alleles)

    assert filter_to_well_normalised(anno_matrix).count_rows() == length
