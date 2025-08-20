"""
unit testing collection for the hail MT methods
"""

import pytest

import hail as hl

from talos.RunHailFiltering import filter_matrix_by_ac, filter_on_quality_flags, populate_callset_frequencies


def test_no_freq_replacement(make_a_mt, caplog):
    """check that when AC/AF/AN are all present, no replacement is attempted"""
    matrix = make_a_mt.annotate_rows(
        info=make_a_mt.info.annotate(
            AC=[1],
            AN=10,
            AF=[0.1],
        ),
    )
    _matrix = populate_callset_frequencies(matrix)
    assert 'AC, AN, AF already present, skipping annotation' in caplog.text, caplog.text


def test_af_replacement(make_a_mt, caplog):
    """check that when AC and AN are present, AF is replaced"""
    matrix = make_a_mt.annotate_rows(
        info=hl.struct(
            AC=[1],
            AN=10,
        ),
    )
    matrix_out = populate_callset_frequencies(matrix)
    assert 'AC, AN present, deriving AF from existing annotations' in caplog.text, caplog.text
    assert matrix_out.info.AF.collect()[0] == [0.1]


def test_full_freq_replacement(make_a_mt, caplog):
    """check that when AC and AN annotations are present, AF is replaced"""
    matrix = make_a_mt.annotate_rows(
        info=hl.struct(),
    )
    matrix_out = populate_callset_frequencies(matrix)
    assert 'Adding AC/AN/AF annotations to MT based on this callset alone' in caplog.text, caplog.text
    assert 'This is unlikely to provide meaningful variant filtering unless this is a huge callset' in caplog.text, (
        caplog.text
    )

    assert matrix_out.info.AF.collect()[0] == [0.5]
    assert matrix_out.info.AC.collect()[0] == [1]
    assert matrix_out.info.AN.collect()[0] == 2


@pytest.mark.parametrize(  # needs clinvar
    'ac,af,clinvar,threshold,rows',
    [
        (1, 0.0, 0, 0.01, 1),
        (6, 0.1, 0, 0.01, 0),
        (6, 0.1, 1, 0.01, 1),
        (50, 0.001, 0, 0.01, 1),
        (50, 0.2, 0, 0.01, 0),
        (50, 0.2, 1, 0.01, 1),
    ],
)
def test_ac_filter_no_filt(
    ac: int,
    af: float,
    clinvar: int,
    threshold: float,
    rows: int,
    make_a_mt: hl.MatrixTable,
):
    """
    run tests on the ac filtering method
    """
    matrix = make_a_mt.annotate_rows(
        info=make_a_mt.info.annotate(
            clinvar_talos=clinvar,
            AC=[ac],
            AF=[af],
        ),
    )

    assert filter_matrix_by_ac(mt=matrix, af_threshold=threshold).count_rows() == rows


@pytest.mark.parametrize(
    'filters,clinvar,length',
    [
        (hl.empty_set(hl.tstr), 0, 1),
        (hl.literal({'fail'}), 0, 0),
        (hl.literal({'fail'}), 1, 1),
        (hl.literal({'VQSR'}), 0, 0),
        (hl.literal({'VQSR'}), 1, 1),
    ],
)
def test_filter_on_quality_flags(
    filters: hl.set,
    clinvar: int,
    length: int,
    make_a_mt: hl.MatrixTable,
):
    """
    annotate filters and run tests
    """
    # to add new alleles, we need to scrub alleles from the key fields
    anno_matrix = make_a_mt.key_rows_by('locus')
    anno_matrix = anno_matrix.annotate_rows(
        filters=filters,
        info=anno_matrix.info.annotate(clinvar_talos=clinvar),
    )
    assert filter_on_quality_flags(anno_matrix).count_rows() == length
