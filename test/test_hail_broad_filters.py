"""
unit testing collection for the hail MT methods
"""


import pytest
import hail as hl

from reanalysis.hail_methods import (
    filter_matrix_by_ac,
    filter_on_quality_flags,
    filter_to_well_normalised,
    filter_to_population_rare,
)


hl_locus = hl.Locus(contig='chr1', position=1, reference_genome='GRCh38')


@pytest.mark.parametrize(
    'ac,an,threshold,rows',
    [
        (1, 1, 0.01, 1),
        (6, 1, 0.01, 0),
        (6, 70, 0.1, 1),
        (50, 999999, 0.01, 1),
        (50, 50, 0.01, 0),
    ],
)
def test_ac_filter_no_filt(ac: int, an: int, threshold: float, rows: int, hail_matrix):
    """
    run tests on the ac filtering method
    """
    matrix = hail_matrix.annotate_rows(info=hail_matrix.info.annotate(AC=ac, AN=an))
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


@pytest.mark.parametrize(
    'exomes,genomes,length',
    [
        (0, 0, 1),
        (1.0, 0, 0),
        (0.04, 0.04, 1),
    ],
)
def test_filter_rows_for_rare(exomes, genomes, length, hail_matrix):
    """
    :param hail_matrix:
    :return:
    """
    anno_matrix = hail_matrix.annotate_rows(
        info=hail_matrix.info.annotate(
            gnomad_ex_af=exomes,
            gnomad_af=genomes,
        )
    )
    matrix = filter_to_population_rare(anno_matrix, thresh=0.05)
    assert matrix.count_rows() == length
