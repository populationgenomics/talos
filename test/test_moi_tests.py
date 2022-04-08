"""
pytests relating to the MOI filters
"""

from typing import List
import pytest
from reanalysis.moi_tests import check_for_second_hit, MOIRunner

from reanalysis.utils import PedPerson


TINY_PEDIGREE = {'test': PedPerson('sample', True, True)}
TINY_CONFIG = {'test': 'test'}
TINY_COMP_HET = {}


@pytest.mark.parametrize(
    'first,comp_hets,sample,gene,truth,values',
    (
        ('', {}, '', '', False, []),  # no values
        ('', {}, 'a', '', False, []),  # sample not present
        ('', {'a': {'c': {'foo': []}}}, 'a', 'b', False, []),  # gene not present
        ('', {'a': {'b': {'foo': []}}}, 'a', 'b', False, []),  # var not present
        (
            'foo',
            {'a': {'b': {'foo': ['bar']}}},
            'a',
            'b',
            True,
            ['bar'],
        ),  # all values present
        (
            'foo',
            {'a': {'b': {'foo': ['bar', 'baz']}}},
            'a',
            'b',
            True,
            ['bar', 'baz'],
        ),  # all values present
    ),
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


@pytest.mark.parametrize(
    'moi_string,filters',
    (
        ('Monoallelic', ['DominantAutosomal']),
        ('Mono_And_Biallelic', ['DominantAutosomal', 'RecessiveAutosomal']),
        ('Unknown', ['DominantAutosomal', 'RecessiveAutosomal']),
        ('Biallelic', ['RecessiveAutosomal']),
        (
            'Hemi_Mono_In_Female',
            ['XRecessive', 'XDominant'],
        ),
        ('Hemi_Bi_In_Female', ['XRecessive']),
        ('Y_Chrom_Variant', ['YHemi']),
    ),
)
def test_moi_runner(moi_string: str, filters: List[str]):
    """

    :param moi_string:
    :param filters:
    :return:
    """
    test_runner = MOIRunner(
        pedigree=TINY_PEDIGREE,
        target_moi=moi_string,
        config=TINY_CONFIG,
        comp_het_lookup=TINY_COMP_HET,
    )
    for filter1, filter2 in zip(test_runner.filter_list, filters):
        assert str(filter1.__class__).__contains__(filter2)
