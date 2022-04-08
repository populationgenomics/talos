"""
tests relating to the MOI filters
"""


from dataclasses import dataclass
from typing import Any, Dict, List, Set
import pytest
from reanalysis.moi_tests import (
    check_for_second_hit,
    BaseMoi,
    DominantAutosomal,
    GNOMAD_AD_AC_THRESHOLD,
    GNOMAD_DOM_HOM_THRESHOLD,
    GNOMAD_RARE_THRESHOLD,
    GNOMAD_REC_HOM_THRESHOLD,
    MOIRunner,
)

from reanalysis.utils import PedPerson


MOI_CONF = {
    GNOMAD_REC_HOM_THRESHOLD: 2,
    GNOMAD_DOM_HOM_THRESHOLD: 1,
    GNOMAD_AD_AC_THRESHOLD: 10,
    GNOMAD_RARE_THRESHOLD: 0.01,
}
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

    # string-comparison
    # the imported (uninstantiated) objects don't have __class__
    # and the instantiated objects don't have a __name__
    for filter1, filter2 in zip(test_runner.filter_list, filters):
        assert str(filter1.__class__).__contains__(filter2)


def test_base_moi():
    """
    test class for the MOI base methods
    :return:
    """

    base = BaseMoi(
        pedigree=TINY_PEDIGREE, config=TINY_CONFIG, applied_moi='test', comp_het={}
    )
    assert base.is_affected('test')


@dataclass
class SimpleVariant:
    """
    a fake version of AbstractVariant
    """

    info: Dict[str, Any]
    het_samples: Set[str]
    hom_samples: Set[str]


def test_dominant_autosomal_passes():
    """
    test case for autosomal dominant
    :return:
    """

    info_dict = {
        'gnomad_af': 0.0001,
        'gnomad_ac': 0,
        'gnomad_hom': 0,
        'gnomad_ex_hom': 0,
        'exac_ac_hom': 0,
    }

    dom = DominantAutosomal(pedigree=TINY_PEDIGREE, config=MOI_CONF)

    # passes with heterozygous
    passing_variant = SimpleVariant(
        info=info_dict, het_samples={'test'}, hom_samples=set()
    )
    results = dom.run(principal_var=passing_variant, ensg='test')
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Dominant'}

    # also passes with homozygous
    passing_variant = SimpleVariant(
        info=info_dict, het_samples=set(), hom_samples={'test'}
    )
    results = dom.run(principal_var=passing_variant, ensg='test')
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Dominant'}

    # no results if no samples
    passing_variant = SimpleVariant(
        info=info_dict, het_samples=set(), hom_samples=set()
    )
    assert not dom.run(principal_var=passing_variant, ensg='test')


@pytest.mark.parametrize(
    'info',
    [
        {'gnomad_af': 0.1},
        {'gnomad_hom': 2},
        {'gnomad_ex_hom': 12},
        {'exac_ac_hom': 12},
    ],
)
def test_dominant_autosomal_fails(info):
    """
    test case for autosomal dominant
    :param info: info dict for the variant
    :return:
    """

    dom = DominantAutosomal(pedigree=TINY_PEDIGREE, config=MOI_CONF)

    # fails due to high af
    failing_variant = SimpleVariant(info=info, het_samples={'test'}, hom_samples=set())
    assert not dom.run(principal_var=failing_variant, ensg='test')
