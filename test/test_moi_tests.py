"""
tests relating to the MOI filters
"""

from unittest import mock

import pytest

from reanalysis.models import Coordinates, SmallVariant
from reanalysis.moi_tests import (
    check_for_second_hit,
    BaseMoi,
    DominantAutosomal,
    MOIRunner,
    RecessiveAutosomalCH,
    RecessiveAutosomalHomo,
    XDominant,
    XRecessiveMale,
    XRecessiveFemaleCH,
    XRecessiveFemaleHom,
)


TEST_COORDS = Coordinates(chrom='1', pos=1, ref='A', alt='C')
TEST_COORDS2 = Coordinates(chrom='2', pos=2, ref='G', alt='T')
TEST_COORDS_X_1 = Coordinates(chrom='X', pos=1, ref='G', alt='T')
TEST_COORDS_X_2 = Coordinates(chrom='X', pos=2, ref='G', alt='T')


@pytest.mark.parametrize(
    'first,comp_hets,sample,values',
    (
        ('', {}, '', []),  # no values
        ('', {}, 'a', []),  # sample not present
        ('', {'a': {'foo': []}}, 'a', []),  # var not present
        (
            'foo',
            {'a': {'foo': ['bar']}},
            'a',
            ['bar'],
        ),  # all values present
        (
            'foo',
            {'a': {'foo': ['bar', 'baz']}},
            'a',
            ['bar', 'baz'],
        ),  # all values present
    ),
)
def test_check_second_hit(first, comp_hets, sample, values):
    """
    quick test for the 2nd hit mechanic
    return all strings when the comp-het lookup contains:
        - the sample
        - the gene
        - the variant signature
    :return:
    """

    assert (
        check_for_second_hit(first_variant=first, comp_hets=comp_hets, sample=sample)
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
    ),
)
def test_moi_runner(moi_string: str, filters: list[str], peddy_ped):
    """

    :param moi_string:
    :param filters:
    :return:
    """
    test_runner = MOIRunner(pedigree=peddy_ped, target_moi=moi_string)

    # string-comparison
    # the imported (uninstantiated) objects don't have __class__
    # and the instantiated objects don't have a __name__
    for filter1, filter2 in zip(test_runner.filter_list, filters):
        assert filter2 in str(filter1.__class__)


def test_dominant_autosomal_fails_on_depth(peddy_ped):
    """
    test case for autosomal dominant depth failure
    :return:
    """

    info_dict = {
        'gnomad_af': 0.0001,
        'gnomad_ac': 0,
        'gnomad_hom': 0,
        'gene_id': 'TEST1',
    }

    dom = DominantAutosomal(pedigree=peddy_ped)

    # passes with heterozygous
    shallow_variant = SmallVariant(
        info=info_dict,
        het_samples={'male'},
        coordinates=TEST_COORDS,
        depths={'male': 1},
        transcript_consequences=[],
    )
    results = dom.run(principal=shallow_variant)  # noqa
    assert len(results) == 0


def test_dominant_autosomal_passes(peddy_ped):
    """
    test case for autosomal dominant
    :return:
    """

    info_dict = {
        'gnomad_af': 0.0001,
        'gnomad_ac': 0,
        'gnomad_hom': 0,
        'cat1': True,
        'gene_id': 'TEST1',
    }

    # attributes relating to categorisation
    boolean_categories = ['cat1']

    dom = DominantAutosomal(pedigree=peddy_ped)

    # passes with heterozygous
    passing_variant = SmallVariant(
        info=info_dict,
        het_samples={'male'},
        coordinates=TEST_COORDS,
        boolean_categories=boolean_categories,
        depths={'male': 999},
        transcript_consequences=[],
    )
    results = dom.run(principal=passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Dominant'}

    # also passes with homozygous
    passing_variant = SmallVariant(
        info=info_dict,
        hom_samples={'male'},
        coordinates=TEST_COORDS,
        boolean_categories=boolean_categories,
        depths={'male': 999},
        transcript_consequences=[],
    )
    results = dom.run(principal=passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Dominant'}

    # no results if no samples
    passing_variant = SmallVariant(
        info=info_dict,
        coordinates=TEST_COORDS,
        boolean_categories=boolean_categories,
        transcript_consequences=[],
    )
    assert len(dom.run(principal=passing_variant)) == 0


@pytest.mark.parametrize(
    'info',
    [
        {'gnomad_af': 0.1},
        {'gnomad_hom': 2},
    ],
)
def test_dominant_autosomal_fails(info, peddy_ped):
    """
    test case for autosomal dominant
    :param info: info dict for the variant
    :return:
    """

    dom = DominantAutosomal(pedigree=peddy_ped)

    # fails due to high af
    failing_variant = SmallVariant(
        info=info,
        het_samples={'male'},
        coordinates=TEST_COORDS,
        transcript_consequences=[],
    )
    assert not dom.run(principal=failing_variant)


def test_recessive_autosomal_hom_passes(peddy_ped):
    """
    check that when the info values are defaults (0)
    we accept a homozygous variant as a Recessive
    """
    passing_variant = SmallVariant(
        hom_samples={'male'},
        coordinates=TEST_COORDS,
        ab_ratios={'male': 1.0},
        depths={'male': 15},
        boolean_categories=['categoryboolean1'],
        info={'categoryboolean1': True, 'gene_id': 'TEST1'},
        transcript_consequences=[],
    )
    rec = RecessiveAutosomalHomo(pedigree=peddy_ped)
    results = rec.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Recessive Homozygous'}


def test_recessive_autosomal_hom_passes_with_ab_flag(peddy_ped):
    """
    check that when the info values are defaults (0)
    we accept a homozygous variant as a Recessive
    """

    passing_variant = SmallVariant(
        hom_samples={'male'},
        coordinates=TEST_COORDS,
        ab_ratios={'male': 0.4},
        depths={'male': 40},
        boolean_categories=['categoryboolean1'],
        info={'categoryboolean1': True, 'gene_id': 'TEST1'},
        transcript_consequences=[],
    )
    rec = RecessiveAutosomalHomo(pedigree=peddy_ped)
    results = rec.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Recessive Homozygous'}
    assert passing_variant.get_sample_flags('male') == {'AB Ratio'}


def test_recessive_autosomal_comp_het_male_passes(peddy_ped):
    """
    check that when the info values are defaults (0)
    and the comp-het test is always True
    we accept a heterozygous variant as a Comp-Het
    """

    passing_variant = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS,
        ab_ratios={'male': 0.5},
        depths={'male': 50},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True},
        transcript_consequences=[],
    )
    passing_variant2 = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS2,
        ab_ratios={'male': 0.5},
        depths={'male': 50},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True},
        transcript_consequences=[],
    )
    comp_hets = {'male': {TEST_COORDS.string_format: [passing_variant2]}}
    rec = RecessiveAutosomalCH(pedigree=peddy_ped)
    results = rec.run(passing_variant, comp_het=comp_hets)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Recessive Comp-Het'}


def test_recessive_autosomal_comp_het_male_passes_partner_flag(peddy_ped):
    """
    check that when the info values are defaults (0)
    and the comp-het test is always True
    we accept a heterozygous variant as a Comp-Het
    """

    passing_variant = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS,
        ab_ratios={'male': 0.5},
        depths={'male': 50},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True},
        transcript_consequences=[],
    )
    passing_variant2 = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS2,
        ab_ratios={'male': 1.0},
        depths={'male': 50},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True},
        transcript_consequences=[],
    )
    comp_hets = {'male': {TEST_COORDS.string_format: [passing_variant2]}}
    rec = RecessiveAutosomalCH(pedigree=peddy_ped)
    results = rec.run(passing_variant, comp_het=comp_hets)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Recessive Comp-Het'}
    assert results[0].flags == {'AB Ratio'}


def test_recessive_autosomal_comp_het_female_passes(peddy_ped):
    """
    check that when the info values are defaults (0)
    and the comp-het test is always True
    we accept a heterozygous variant as a Comp-Het
    :return:
    """

    passing_variant = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS,
        ab_ratios={'female': 0.5},
        depths={'female': 50},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True},
        transcript_consequences=[],
    )
    passing_variant2 = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS2,
        ab_ratios={'female': 0.5},
        depths={'female': 50},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True},
        transcript_consequences=[],
    )
    comp_hets = {'female': {TEST_COORDS.string_format: [passing_variant2]}}
    rec = RecessiveAutosomalCH(pedigree=peddy_ped)
    results = rec.run(passing_variant, comp_het=comp_hets)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Recessive Comp-Het'}
    assert not results[0].flags


def test_recessive_autosomal_comp_het_fails_no_ch_return(peddy_ped):
    """
    check that when the info values are defaults (0)
    and the comp-het test is always False
    we have no accepted MOI

    :return:
    """

    failing_variant = SmallVariant(
        info={'gene_id': 'TEST1'},
        het_samples={'male'},
        depths={'male': 50},
        coordinates=TEST_COORDS,
        transcript_consequences=[],
    )
    rec = RecessiveAutosomalCH(pedigree=peddy_ped)
    assert not rec.run(failing_variant)


def test_recessive_autosomal_comp_het_fails_no_paired_call(peddy_ped):
    """
    check that when the info values are defaults (0)
    and the comp-het test is always False
    we have no accepted MOI

    :return:
    """

    failing_variant = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS,
        ab_ratios={'male': 0.5},
        depths={'male': 50},
        info={'gene_id': 'TEST1'},
        transcript_consequences=[],
    )
    failing_variant2 = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS2,
        ab_ratios={'female': 0.5},
        depths={'female': 50},
        info={'gene_id': 'TEST1'},
        transcript_consequences=[],
    )

    rec = RecessiveAutosomalCH(pedigree=peddy_ped)
    assert not rec.run(
        failing_variant,
        comp_het={'male': {TEST_COORDS2.string_format: [failing_variant2]}},
    )


@pytest.mark.parametrize(
    'info', [{'gnomad_hom': 3, 'gene_id': 'TEST1'}]
)  # threshold is 2
def test_recessive_autosomal_hom_fails(info, peddy_ped):
    """
    check that when the info values are failures
    we have no confirmed MOI
    """

    failing_variant = SmallVariant(
        info=info,
        het_samples={'male'},
        hom_samples={'male'},
        coordinates=TEST_COORDS,
        transcript_consequences=[],
    )
    rec = RecessiveAutosomalHomo(pedigree=peddy_ped)
    assert not rec.run(failing_variant)


def test_x_dominant_female_and_male_het_passes(peddy_ped):
    """
    check that a male is accepted as a het
    :return:
    """
    passing_variant = SmallVariant(
        boolean_categories=['categoryboolean1'],
        info={'gnomad_hemi': 0, 'gene_id': 'TEST1', 'categoryboolean1': True},
        het_samples={'female', 'male'},
        depths={'female': 50, 'male': 50},
        coordinates=TEST_COORDS_X_1,
        transcript_consequences=[],
    )
    x_dom = XDominant(pedigree=peddy_ped)
    results = x_dom.run(passing_variant)

    assert len(results) == 2
    reasons = {result.reasons.pop() for result in results}
    assert reasons == {'X_Dominant', 'X_Dominant'}


def test_x_dominant_female_hom_passes(peddy_ped):
    """
    check that a male is accepted as a het
    :return:
    """
    passing_variant = SmallVariant(
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True},
        hom_samples={'female'},
        depths={'female': 100},
        ab_ratios={'female': 0.5},
        coordinates=TEST_COORDS_X_1,
        transcript_consequences=[],
    )
    x_dom = XDominant(pedigree=peddy_ped)
    results = x_dom.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'X_Dominant'}


def test_x_dominant_male_hom_passes(peddy_ped):
    """
    check that a male is accepted as a het
    :return:
    """
    passing_variant = SmallVariant(
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True},
        hom_samples={'male'},
        depths={'male': 100},
        coordinates=TEST_COORDS_X_1,
        transcript_consequences=[],
    )
    x_dom = XDominant(pedigree=peddy_ped)
    results = x_dom.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'X_Dominant'}


@pytest.mark.parametrize(
    'info',
    [
        {'gnomad_af': 0.1, 'gene_id': 'TEST1', 'categoryboolean1': True},
        {'gnomad_hom': 2, 'gene_id': 'TEST1', 'categoryboolean1': True},
        {'gnomad_hemi': 3, 'gene_id': 'TEST1', 'categoryboolean1': True},
    ],
)
def test_x_dominant_info_fails(info, peddy_ped):
    """
    check for info dict exclusions
    """
    passing_variant = SmallVariant(
        info=info,
        hom_samples={'male'},
        coordinates=TEST_COORDS_X_1,
        boolean_categories=['categoryboolean1'],
        transcript_consequences=[],
    )
    x_dom = XDominant(pedigree=peddy_ped)
    assert len(x_dom.run(passing_variant)) == 0


def test_x_recessive_male_hom_passes(peddy_ped):
    """

    :return:
    """

    passing_variant = SmallVariant(
        hom_samples={'female', 'male'},
        coordinates=TEST_COORDS_X_1,
        ab_ratios={'female': 1.0, 'male': 1.0},
        depths={'female': 100, 'male': 100},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True},
        transcript_consequences=[],
    )
    x_rec = XRecessiveMale(pedigree=peddy_ped)
    results = x_rec.run(passing_variant, comp_het={})
    assert len(results) == 1
    assert results[0].reasons == {'X_Male'}


def test_x_recessive_female_hom_passes(peddy_ped):
    """
    :return:
    """

    passing_variant = SmallVariant(
        hom_samples={'female', 'male'},
        coordinates=TEST_COORDS_X_1,
        ab_ratios={'female': 1.0, 'male': 1.0},
        depths={'female': 100, 'male': 100},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True},
        transcript_consequences=[],
    )
    x_rec = XRecessiveFemaleHom(pedigree=peddy_ped)
    results = x_rec.run(passing_variant, comp_het={})
    assert len(results) == 1
    assert results[0].reasons == {'X_Recessive HOM Female'}


def test_x_recessive_male_het_passes(peddy_ped):
    """

    :return:
    """
    passing_variant = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS_X_1,
        ab_ratios={'male': 0.5},
        depths={'male': 50},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True},
        transcript_consequences=[],
    )
    x_rec = XRecessiveMale(pedigree=peddy_ped)
    results = x_rec.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'X_Male'}


def test_x_recessive_female_het_passes(peddy_ped):
    """

    :return:
    """

    passing_variant = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_1,
        ab_ratios={'female': 0.5},
        depths={'female': 50},
        sample_categories=['categorysample4'],
        info={
            'gene_id': 'TEST1',
            'categorysample4': ['female'],
        },
        transcript_consequences=[],
    )
    passing_variant_2 = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_2,
        ab_ratios={'female': 0.5},
        depths={'female': 50},
        sample_categories=['categorysample4'],
        info={
            'gene_id': 'TEST1',
            'categorysample4': ['female'],
        },
        transcript_consequences=[],
    )
    comp_hets = {'female': {'X-1-G-T': [passing_variant_2]}}
    x_rec = XRecessiveFemaleCH(pedigree=peddy_ped)
    results = x_rec.run(passing_variant, comp_het=comp_hets)
    assert len(results) == 1
    assert results[0].reasons == {'X_RecessiveFemaleCompHet'}


def test_het_de_novo_passes(peddy_ped):
    """

    :return:
    """

    passing_variant = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_1,
        sample_categories=['categorysample4'],
        ab_ratios={'female': 0.5},
        depths={'female': 99},
        info={'gene_id': 'TEST1', 'categorysample4': ['female']},
        transcript_consequences=[],
    )
    dom_a = DominantAutosomal(pedigree=peddy_ped)
    results = dom_a.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Dominant'}
    assert not results[0].flags


def test_het_de_novo_het_passes_flagged(peddy_ped):
    """

    :return:
    """

    passing_variant = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_1,
        sample_categories=['categorysample4'],
        ab_ratios={'female': 0.5},
        depths={'female': 99},
        info={'gene_id': 'TEST1', 'categorysample4': ['female']},
        transcript_consequences=[],
    )
    dom_a = DominantAutosomal(pedigree=peddy_ped)
    results = dom_a.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Dominant'}


def test_x_recessive_female_het_fails(peddy_ped):
    """
    :return:
    """

    passing_variant = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_1,
        ab_ratios={'female': 0.5},
        depths={'female': 50},
        sample_categories=['categorysample4'],
        info={
            'gene_id': 'TEST1',
            'categorysample4': ['male'],
        },
        transcript_consequences=[],
    )
    passing_variant_2 = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS_X_2,
        ab_ratios={'male': 0.5},
        depths={'male': 50},
        sample_categories=['categorysample4'],
        info={
            'gene_id': 'TEST1',
            'categorysample4': ['male'],
        },
        transcript_consequences=[],
    )
    comp_hets = {'female': {'x-2-A-C': [passing_variant_2]}}
    x_rec = XRecessiveFemaleCH(pedigree=peddy_ped)
    results = x_rec.run(passing_variant, comp_het=comp_hets)
    assert not results


@mock.patch('reanalysis.moi_tests.check_for_second_hit')
def test_x_recessive_female_het_no_pair_fails(second_hit: mock.Mock, peddy_ped):
    """ """

    passing_variant = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_1,
        ab_ratios={'female': 0.5},
        depths={'female': 50},
        info={
            'gene_id': 'TEST1',
            'categorysample1': True,
            'boolean_categories': 'categorysample1',
        },
        transcript_consequences=[],
    )
    second_hit.return_value = []
    x_rec = XRecessiveFemaleCH(pedigree=peddy_ped)
    assert not x_rec.run(passing_variant)


# trio male, mother_1, father_1; only 'male' is affected
def test_check_familial_inheritance_simple(peddy_ped):
    """
    test the check_familial_inheritance method
    :return:
    """

    base_moi = BaseMoi(pedigree=peddy_ped, applied_moi='applied')

    result = base_moi.check_familial_inheritance(
        sample_id='male', called_variants={'male'}
    )
    assert result


def test_check_familial_inheritance_mother_fail(peddy_ped):
    """
    test the check_familial_inheritance method
    :return:
    """

    base_moi = BaseMoi(pedigree=peddy_ped, applied_moi='applied')

    result = base_moi.check_familial_inheritance(
        sample_id='male', called_variants={'male', 'mother_1'}
    )
    assert not result


def test_check_familial_inheritance_mother_passes(peddy_ped):
    """
    test the check_familial_inheritance method
    mother in variant calls, but partial penetrance
    :return:
    """

    base_moi = BaseMoi(pedigree=peddy_ped, applied_moi='applied')

    result = base_moi.check_familial_inheritance(
        sample_id='male',
        called_variants={'male', 'mother_1'},
        partial_pen=True,
    )
    assert result


def test_check_familial_inheritance_father_fail(peddy_ped):
    """
    test the check_familial_inheritance method
    :return:
    """

    base_moi = BaseMoi(pedigree=peddy_ped, applied_moi='applied')

    result = base_moi.check_familial_inheritance(
        sample_id='male', called_variants={'male', 'father_1'}
    )
    assert not result


def test_check_familial_inheritance_father_passes(peddy_ped):
    """
    test the check_familial_inheritance method
    father in variant calls, but partial penetrance
    :return:
    """

    base_moi = BaseMoi(pedigree=peddy_ped, applied_moi='applied')

    result = base_moi.check_familial_inheritance(
        sample_id='male',
        called_variants={'male', 'father_1'},
        partial_pen=True,
    )
    assert result


def test_check_familial_inheritance_top_down(peddy_ped):
    """
    test the check_familial_inheritance method
    father in variant calls, but partial penetrance
    :return:
    """

    base_moi = BaseMoi(pedigree=peddy_ped, applied_moi='applied')

    result = base_moi.check_familial_inheritance(
        sample_id='father_1',
        called_variants={'male', 'father_1'},
        partial_pen=True,
    )
    assert result


def test_check_familial_inheritance_no_calls(peddy_ped):
    """
    test the check_familial_inheritance method where there are no calls
    will fail as affected proband not in calls
    """

    base_moi = BaseMoi(pedigree=peddy_ped, applied_moi='applied')

    result = base_moi.check_familial_inheritance(
        sample_id='male',
        called_variants=set(),
        partial_pen=True,
    )
    assert not result


def test_genotype_calls(peddy_ped):
    """
    test the manual genotype assignments
    """
    base_moi = DominantAutosomal(pedigree=peddy_ped, applied_moi='applied')

    info_dict = {
        'gnomad_af': 0.0001,
        'gnomad_ac': 0,
        'gnomad_hom': 0,
        'gene_id': 'TEST1',
    }
    variant = SmallVariant(
        info=info_dict,
        het_samples={'male'},
        hom_samples={'female'},
        coordinates=TEST_COORDS,
        transcript_consequences=[],
    )
    assert base_moi.get_family_genotypes(variant, 'male') == {
        'father_1': 'WT',
        'male': 'Het',
        'mother_1': 'WT',
    }
    assert base_moi.get_family_genotypes(variant, 'female') == {
        'father_2': 'WT',
        'female': 'Hom',
        'mother_2': 'WT',
    }
    x_variant = SmallVariant(
        info=info_dict,
        het_samples={'male', 'female'},
        coordinates=TEST_COORDS_X_1,
        transcript_consequences=[],
    )
    assert base_moi.get_family_genotypes(x_variant, 'male') == {
        'father_1': 'WT',
        'male': 'Hemi',
        'mother_1': 'WT',
    }
    assert base_moi.get_family_genotypes(x_variant, 'female') == {
        'father_2': 'WT',
        'female': 'Het',
        'mother_2': 'WT',
    }

    x_variant_2 = SmallVariant(
        info=info_dict,
        hom_samples={'male', 'female'},
        coordinates=TEST_COORDS_X_1,
        transcript_consequences=[],
    )
    assert base_moi.get_family_genotypes(x_variant_2, 'male') == {
        'father_1': 'WT',
        'male': 'Hemi',
        'mother_1': 'WT',
    }
    assert base_moi.get_family_genotypes(x_variant_2, 'female') == {
        'father_2': 'WT',
        'female': 'Hom',
        'mother_2': 'WT',
    }

    variant_missing = SmallVariant(
        info=info_dict,
        coordinates=TEST_COORDS,
        transcript_consequences=[],
    )
    assert base_moi.get_family_genotypes(variant_missing, 'male') == {
        'father_1': 'WT',
        'male': 'WT',
        'mother_1': 'WT',
    }
