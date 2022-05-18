"""
tests relating to the MOI filters
"""


from dataclasses import dataclass, field
from typing import Any, Dict, List, Set

import os

from unittest import mock
from peddy.peddy import Ped

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
    RecessiveAutosomal,
    XDominant,
    XRecessive,
)

from reanalysis.utils import Coordinates


PWD = os.path.dirname(__file__)
INPUT = os.path.join(PWD, 'input')
PED_FILE = os.path.join(INPUT, 'pedfile.ped')

MOI_CONF = {
    GNOMAD_REC_HOM_THRESHOLD: 2,
    GNOMAD_DOM_HOM_THRESHOLD: 1,
    GNOMAD_AD_AC_THRESHOLD: 10,
    GNOMAD_RARE_THRESHOLD: 0.01,
}
TEST_COORDS = Coordinates('1', 1, 'A', 'C')
TEST_COORDS2 = Coordinates('2', 2, 'G', 'T')
TINY_PEDIGREE = Ped(PED_FILE)
TINY_CONFIG = {'male': 'male'}
TINY_COMP_HET = {}


@dataclass
class SimpleVariant:
    """
    a fake version of AbstractVariant
    """

    info: Dict[str, Any]
    het_samples: Set[str]
    hom_samples: Set[str]
    coords: Coordinates
    category_1: bool = True
    category_4: List[str] = field(default_factory=list)

    def sample_specific_category_check(self, sample):
        """
        pass
        :param sample:
        :return:
        """
        return sample in self.category_4


@dataclass
class RecessiveSimpleVariant:
    """
    a fake version of AbstractVariant
    """

    info: Dict[str, Any]
    het_samples: Set[str]
    hom_samples: Set[str]
    coords: Coordinates
    category_4: List[str]
    # add category default
    category_1: bool = True

    @property
    def category_1_2_3(self):
        """
        mock method
        :return:
        """
        return self.category_1

    def sample_de_novo(self, sample):
        """
        pass
        :param sample:
        :return:
        """
        return sample in self.category_4

    def sample_specific_category_check(self, sample):
        """
        pass
        :param sample:
        :return:
        """
        return (sample in self.category_4) or self.category_1_2_3


@pytest.mark.parametrize(
    'first,comp_hets,sample,gene,values',
    (
        ('', {}, '', '', []),  # no values
        ('', {}, 'a', '', []),  # sample not present
        ('', {'a': {'c': {'foo': []}}}, 'a', 'b', []),  # gene not present
        ('', {'a': {'b': {'foo': []}}}, 'a', 'b', []),  # var not present
        (
            'foo',
            {'a': {'b': {'foo': ['bar']}}},
            'a',
            'b',
            ['bar'],
        ),  # all values present
        (
            'foo',
            {'a': {'b': {'foo': ['bar', 'baz']}}},
            'a',
            'b',
            ['bar', 'baz'],
        ),  # all values present
    ),
)
def test_check_second_hit(first, comp_hets, sample, gene, values):
    """
    quick test for the 2nd hit mechanic
    return all strings when the comp-het lookup contains:
        - the sample
        - the gene
        - the variant signature
    :return:
    """

    assert (
        check_for_second_hit(
            first_variant=first, comp_hets=comp_hets, sample=sample, gene=gene
        )
        == values
    )


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
        info=info_dict, het_samples={'male'}, hom_samples=set(), coords=TEST_COORDS
    )
    results = dom.run(principal_var=passing_variant, gene_lookup={})
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Dominant'}

    # also passes with homozygous
    passing_variant = SimpleVariant(
        info=info_dict, het_samples=set(), hom_samples={'male'}, coords=TEST_COORDS
    )
    results = dom.run(principal_var=passing_variant, gene_lookup={})
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Dominant'}

    # no results if no samples
    passing_variant = SimpleVariant(
        info=info_dict, het_samples=set(), hom_samples=set(), coords=TEST_COORDS
    )
    assert not dom.run(principal_var=passing_variant, gene_lookup={})


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
    failing_variant = SimpleVariant(
        info=info, het_samples={'male'}, hom_samples=set(), coords=TEST_COORDS
    )
    assert not dom.run(principal_var=failing_variant, gene_lookup={})


def test_recessive_autosomal_hom_passes():
    """
    check that when the info values are defaults (0)
    we accept a homozygous variant as a Recessive
    """

    passing_variant = RecessiveSimpleVariant(
        info={},
        het_samples=set(),
        hom_samples={'male'},
        coords=TEST_COORDS,
        category_4=[],
    )
    rec = RecessiveAutosomal(
        pedigree=TINY_PEDIGREE, config={GNOMAD_REC_HOM_THRESHOLD: 1}, comp_het={}
    )
    results = rec.run(passing_variant, gene_lookup={})
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Recessive Homozygous'}


@mock.patch('reanalysis.moi_tests.check_for_second_hit')
def test_recessive_autosomal_comp_het_male_passes(second_hit: mock.Mock):
    """
    check that when the info values are defaults (0)
    and the comp-het test is always True
    we accept a heterozygous variant as a Comp-Het

    :param second_hit: patch
    :return:
    """

    passing_variant = RecessiveSimpleVariant(
        info={},
        het_samples={'male'},
        hom_samples=set(),
        coords=TEST_COORDS,
        category_4=[],
    )
    passing_variant2 = RecessiveSimpleVariant(
        info={},
        het_samples={'male'},
        hom_samples=set(),
        coords=TEST_COORDS2,
        category_4=[],
    )
    second_hit.return_value = [TEST_COORDS2.string_format]
    gene_var_dict = {TEST_COORDS2.string_format: passing_variant2}
    rec = RecessiveAutosomal(
        pedigree=TINY_PEDIGREE, config={GNOMAD_REC_HOM_THRESHOLD: 1}, comp_het={}
    )
    results = rec.run(passing_variant, gene_lookup=gene_var_dict)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Recessive Compound-Het'}


@mock.patch('reanalysis.moi_tests.check_for_second_hit')
def test_recessive_autosomal_comp_het_female_passes(second_hit: mock.Mock):
    """
    check that when the info values are defaults (0)
    and the comp-het test is always True
    we accept a heterozygous variant as a Comp-Het

    :param second_hit: patch
    :return:
    """

    passing_variant = RecessiveSimpleVariant(
        info={},
        het_samples={'female'},
        hom_samples=set(),
        coords=TEST_COORDS,
        category_4=[],
    )
    passing_variant2 = RecessiveSimpleVariant(
        info={},
        het_samples={'female'},
        hom_samples=set(),
        coords=TEST_COORDS2,
        category_4=[],
    )
    second_hit.return_value = [TEST_COORDS2.string_format]
    gene_var_dict = {TEST_COORDS2.string_format: passing_variant2}
    rec = RecessiveAutosomal(
        pedigree=TINY_PEDIGREE, config={GNOMAD_REC_HOM_THRESHOLD: 1}, comp_het={}
    )
    results = rec.run(passing_variant, gene_lookup=gene_var_dict)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Recessive Compound-Het'}


@mock.patch('reanalysis.moi_tests.check_for_second_hit')
def test_recessive_autosomal_comp_het_fails_no_ch_return(second_hit: mock.Mock):
    """
    check that when the info values are defaults (0)
    and the comp-het test is always False
    we have no accepted MOI

    :param second_hit: patch
    :return:
    """

    failing_variant = SimpleVariant(
        info={}, het_samples={'male'}, hom_samples=set(), coords=TEST_COORDS
    )
    second_hit.return_value = []
    rec = RecessiveAutosomal(
        pedigree=TINY_PEDIGREE, config={GNOMAD_REC_HOM_THRESHOLD: 1}, comp_het={}
    )
    assert not rec.run(failing_variant, gene_lookup={})


@mock.patch('reanalysis.moi_tests.check_for_second_hit')
def test_recessive_autosomal_comp_het_fails_no_paired_call(second_hit: mock.Mock):
    """
    check that when the info values are defaults (0)
    and the comp-het test is always False
    we have no accepted MOI

    :param second_hit: patch
    :return:
    """

    failing_variant = RecessiveSimpleVariant(
        info={},
        het_samples={'male'},
        hom_samples=set(),
        coords=TEST_COORDS,
        category_4=[],
    )
    failing_variant2 = RecessiveSimpleVariant(
        info={},
        het_samples={'female'},
        hom_samples=set(),
        coords=TEST_COORDS2,
        category_4=[],
    )
    second_hit.return_value = [TEST_COORDS2.string_format]
    gene_var_dict = {TEST_COORDS2.string_format: failing_variant2}

    rec = RecessiveAutosomal(
        pedigree=TINY_PEDIGREE, config={GNOMAD_REC_HOM_THRESHOLD: 1}, comp_het={}
    )
    assert not rec.run(failing_variant, gene_lookup=gene_var_dict)


@pytest.mark.parametrize(
    'info',
    [
        {'gnomad_hom': 2},
        {'gnomad_ex_hom': 12},
        {'exac_ac_hom': 12},
    ],
)
def test_recessive_autosomal_hom_fails(info):
    """
    check that when the info values are failures
    we have no confirmed MOI
    """

    failing_variant = SimpleVariant(
        info=info, het_samples={'male'}, hom_samples={'male'}, coords=TEST_COORDS
    )
    rec = RecessiveAutosomal(
        pedigree=TINY_PEDIGREE, config={GNOMAD_REC_HOM_THRESHOLD: 1}, comp_het={}
    )
    assert not rec.run(failing_variant, gene_lookup={})


def test_x_dominant_female_and_male_het_passes():
    """
    check that a male is accepted as a het
    :return:
    """
    x_coords = Coordinates('x', 1, 'A', 'C')
    passing_variant = SimpleVariant(
        info={}, het_samples={'female', 'male'}, hom_samples=set(), coords=x_coords
    )
    x_dom = XDominant(pedigree=TINY_PEDIGREE, config=MOI_CONF)
    results = x_dom.run(passing_variant, gene_lookup={})

    assert len(results) == 2
    reasons = sorted([result.reasons.pop() for result in results])
    assert reasons == ['X_Dominant Female', 'X_Dominant Male']


def test_x_dominant_female_hom_passes():
    """
    check that a male is accepted as a het
    :return:
    """
    x_coords = Coordinates('x', 1, 'A', 'C')
    passing_variant = SimpleVariant(
        info={}, hom_samples={'female'}, het_samples=set(), coords=x_coords
    )
    x_dom = XDominant(pedigree=TINY_PEDIGREE, config=MOI_CONF)
    results = x_dom.run(passing_variant, gene_lookup={})
    assert len(results) == 1
    assert results[0].reasons == {'X_Dominant Female'}


def test_x_dominant_male_hom_passes():
    """
    check that a male is accepted as a het
    :return:
    """
    x_coords = Coordinates('x', 1, 'A', 'C')
    passing_variant = SimpleVariant(
        info={}, hom_samples={'male'}, het_samples=set(), coords=x_coords
    )
    x_dom = XDominant(pedigree=TINY_PEDIGREE, config=MOI_CONF)
    results = x_dom.run(passing_variant, gene_lookup={})
    assert len(results) == 1
    assert results[0].reasons == {'X_Dominant Male'}


def test_x_moi_on_non_x_fails():
    """
    check that a male is accepted as a het
    :return:
    """
    y_coords = Coordinates('y', 1, 'A', 'C')
    y_variant = SimpleVariant(
        info={}, het_samples={'male'}, hom_samples=set(), coords=y_coords
    )
    x_dom = XDominant(pedigree=TINY_PEDIGREE, config=MOI_CONF)
    with pytest.raises(Exception):
        x_dom.run(y_variant, gene_lookup={})


@pytest.mark.parametrize(
    'info',
    [
        {'gnomad_af': 0.1},
        {'gnomad_hom': 2},
        {'gnomad_ex_hom': 12},
        {'exac_ac_hom': 12},
    ],
)
def test_x_dominant_info_fails(info):
    """
    check for info dict exclusions
    :param info:
    :return:
    """
    x_coords = Coordinates('x', 1, 'A', 'C')
    passing_variant = SimpleVariant(
        info=info, hom_samples={'female'}, het_samples=set(), coords=x_coords
    )
    x_dom = XDominant(pedigree=TINY_PEDIGREE, config=MOI_CONF)
    assert not x_dom.run(passing_variant, gene_lookup={})


def test_x_recessive_male_and_female_hom_passes():
    """

    :return:
    """

    x_coords = Coordinates('x', 1, 'A', 'C')
    passing_variant = RecessiveSimpleVariant(
        info={},
        hom_samples={'female', 'male'},
        het_samples=set(),
        coords=x_coords,
        category_4=[],
    )
    x_rec = XRecessive(pedigree=TINY_PEDIGREE, config=MOI_CONF, comp_het={})
    results = x_rec.run(passing_variant, gene_lookup={})
    assert len(results) == 2

    reasons = sorted([result.reasons.pop() for result in results])
    assert reasons == ['X_Recessive Female', 'X_Recessive Male']


def test_x_recessive_male_het_passes():
    """

    :return:
    """
    x_coords = Coordinates('x', 1, 'A', 'C')
    passing_variant = RecessiveSimpleVariant(
        info={}, het_samples={'male'}, hom_samples=set(), coords=x_coords, category_4=[]
    )
    x_rec = XRecessive(pedigree=TINY_PEDIGREE, config=MOI_CONF, comp_het={})
    results = x_rec.run(passing_variant, gene_lookup={})
    assert len(results) == 1
    assert results[0].reasons == {'X_Recessive Male'}


def test_x_recessive_y_variant_fails():
    """

    :return:
    """
    y_coords = Coordinates('y', 1, 'A', 'C')
    passing_variant = RecessiveSimpleVariant(
        info={}, hom_samples={'male'}, het_samples=set(), coords=y_coords, category_4=[]
    )
    x_rec = XRecessive(pedigree=TINY_PEDIGREE, config=MOI_CONF, comp_het={})
    with pytest.raises(Exception):
        x_rec.run(passing_variant, gene_lookup={})


@mock.patch('reanalysis.moi_tests.check_for_second_hit')
def test_x_recessive_female_het_passes(second_hit: mock.patch):
    """

    :return:
    """

    second_hit.return_value = ['x-2-A-C']
    passing_variant = RecessiveSimpleVariant(
        info={},
        het_samples={'female'},
        hom_samples=set(),
        coords=Coordinates('x', 1, 'A', 'C'),
        category_4=['female'],
    )
    passing_variant_2 = RecessiveSimpleVariant(
        info={},
        het_samples={'female'},
        hom_samples=set(),
        coords=Coordinates('x', 2, 'A', 'C'),
        category_4=['female'],
    )
    gene_dict = {'x-2-A-C': passing_variant_2}
    x_rec = XRecessive(pedigree=TINY_PEDIGREE, config=MOI_CONF, comp_het={})
    results = x_rec.run(passing_variant, gene_lookup=gene_dict)
    assert len(results) == 1
    assert results[0].reasons == {'X_Recessive Compound-Het Female'}


@mock.patch('reanalysis.moi_tests.check_for_second_hit')
def test_x_recessive_female_het_fails(second_hit: mock.patch):
    """

    :return:
    """

    second_hit.return_value = ['x-2-A-C']
    passing_variant = RecessiveSimpleVariant(
        info={},
        het_samples={'female'},
        hom_samples=set(),
        coords=Coordinates('x', 1, 'A', 'C'),
        category_4=['male'],
    )
    passing_variant_2 = RecessiveSimpleVariant(
        info={},
        het_samples={'male'},
        hom_samples=set(),
        coords=Coordinates('x', 2, 'A', 'C'),
        category_4=['male'],
    )
    gene_dict = {'x-2-A-C': passing_variant_2}
    x_rec = XRecessive(pedigree=TINY_PEDIGREE, config=MOI_CONF, comp_het={})
    results = x_rec.run(passing_variant, gene_lookup=gene_dict)
    assert not results


@mock.patch('reanalysis.moi_tests.check_for_second_hit')
def test_x_recessive_female_het_no_pair_fails(second_hit: mock.patch):
    """

    :return:
    """

    second_hit.return_value = []
    passing_variant = RecessiveSimpleVariant(
        info={},
        het_samples={'female'},
        hom_samples=set(),
        coords=Coordinates('x', 1, 'A', 'C'),
        category_4=[],
    )
    x_rec = XRecessive(pedigree=TINY_PEDIGREE, config=MOI_CONF, comp_het={})
    assert not x_rec.run(passing_variant, gene_lookup={})


# trio male, mother_1, father_1; only 'male' is affected
def test_check_familial_inheritance_simple():
    """
    test the check_familial_inheritance method
    :return:
    """

    base_moi = BaseMoi(
        pedigree=TINY_PEDIGREE, config=TINY_CONFIG, applied_moi='applied', comp_het={}
    )

    result = base_moi.check_familial_inheritance(
        sample_id='male', called_variants={'male'}
    )
    assert result


def test_check_familial_inheritance_mother_fail():
    """
    test the check_familial_inheritance method
    :return:
    """

    base_moi = BaseMoi(
        pedigree=TINY_PEDIGREE, config=TINY_CONFIG, applied_moi='applied', comp_het={}
    )

    result = base_moi.check_familial_inheritance(
        sample_id='male', called_variants={'male', 'mother_1'}
    )
    assert not result


def test_check_familial_inheritance_mother_passes():
    """
    test the check_familial_inheritance method
    mother in variant calls, but partial penetrance
    :return:
    """

    base_moi = BaseMoi(
        pedigree=TINY_PEDIGREE, config=TINY_CONFIG, applied_moi='applied', comp_het={}
    )

    result = base_moi.check_familial_inheritance(
        sample_id='male',
        called_variants={'male', 'mother_1'},
        complete_penetrance=False,
    )
    assert result


def test_check_familial_inheritance_father_fail():
    """
    test the check_familial_inheritance method
    :return:
    """

    base_moi = BaseMoi(
        pedigree=TINY_PEDIGREE, config=TINY_CONFIG, applied_moi='applied', comp_het={}
    )

    result = base_moi.check_familial_inheritance(
        sample_id='male', called_variants={'male', 'father_1'}
    )
    assert not result


def test_check_familial_inheritance_father_passes():
    """
    test the check_familial_inheritance method
    father in variant calls, but partial penetrance
    :return:
    """

    base_moi = BaseMoi(
        pedigree=TINY_PEDIGREE, config=TINY_CONFIG, applied_moi='applied', comp_het={}
    )

    result = base_moi.check_familial_inheritance(
        sample_id='male',
        called_variants={'male', 'father_1'},
        complete_penetrance=False,
    )
    assert result


def test_check_familial_inheritance_top_down():
    """
    test the check_familial_inheritance method
    father in variant calls, but partial penetrance
    :return:
    """

    base_moi = BaseMoi(
        pedigree=TINY_PEDIGREE, config=TINY_CONFIG, applied_moi='applied', comp_het={}
    )

    result = base_moi.check_familial_inheritance(
        sample_id='father_1',
        called_variants={'male', 'father_1'},
        complete_penetrance=False,
    )
    assert result


def test_check_familial_inheritance_no_calls():
    """
    test the check_familial_inheritance method where there are no calls
    will fail as affected proband not in calls
    :return:
    """

    base_moi = BaseMoi(
        pedigree=TINY_PEDIGREE, config=TINY_CONFIG, applied_moi='applied', comp_het={}
    )

    result = base_moi.check_familial_inheritance(
        sample_id='male',
        called_variants=set(),
        complete_penetrance=False,
    )
    # should fail immediately
    assert not result
