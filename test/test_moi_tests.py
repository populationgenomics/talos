"""
pytests relating to the MOI filters
"""

import pytest
from reanalysis.moi_tests import check_for_second_hit


@pytest.mark.parametrize(
    'first,comp_hets,sample,gene,truth,values',
    [
        ['', {}, '', '', False, []],  # no values
        ['', {}, 'a', '', False, []],  # sample not present
        ['', {'a': {'c': {'foo': []}}}, 'a', 'b', False, []],  # gene not present
        ['', {'a': {'b': {'foo': []}}}, 'a', 'b', False, []],  # var not present
        [
            'foo',
            {'a': {'b': {'foo': ['bar']}}},
            'a',
            'b',
            True,
            ['bar'],
        ],  # all values present
        [
            'foo',
            {'a': {'b': {'foo': ['bar', 'baz']}}},
            'a',
            'b',
            True,
            ['bar', 'baz'],
        ],  # all values present
    ],
)
def test_check_second_hit(first, comp_hets, sample, gene, truth, values):
    """
    quick test for the 2nd hit mechanic
    return all strings when the comp-het lookup contains:
        - the sample
        - the gene
        - the variant signature
    :return:
    """

    assert check_for_second_hit(
        first_variant=first, comp_hets=comp_hets, sample=sample, gene=gene
    ) == (truth, values)
