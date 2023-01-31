"""
unit testing collection for the hail MT methods
"""


import pytest
import hail as hl

from reanalysis.hail_filter_and_label import (
    filter_matrix_by_ac,
    filter_on_quality_flags,
    filter_to_well_normalised,
)


hl_locus = hl.Locus(contig='chr1', position=1, reference_genome='GRCh38')


@pytest.mark.parametrize(  # needs clinvar
    'ac,an,clinvar,threshold,rows',
    [
        (1, 1, 0, 0.01, 1),
        (6, 1, 0, 0.01, 0),
        (6, 1, 1, 0.01, 1),
        (6, 70, 0, 0.1, 1),
        (50, 999999, 0, 0.01, 1),
        (50, 50, 0, 0.01, 0),
        (50, 50, 1, 0.01, 1),
    ],
)
def test_ac_filter_no_filt(
    ac: int, an: int, clinvar: int, threshold: float, rows: int, hail_matrix
):
    """
    run tests on the ac filtering method
    check that a clinvar pathogenic overrides the AC test
    """
    matrix = hail_matrix.annotate_rows(AC=ac, AN=an)
    matrix = matrix.annotate_rows(info=matrix.info.annotate(clinvar_aip=clinvar))

    assert filter_matrix_by_ac(matrix, threshold).count_rows() == rows


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
