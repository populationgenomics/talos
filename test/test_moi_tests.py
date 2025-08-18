"""
tests relating to the MOI filters
"""

import pytest

from talos.models import Coordinates, SmallVariant
from talos.moi_tests import (
    BaseMoi,
    ClinVarFilter,
    DominantAutosomal,
    MOIRunner,
    RecessiveAutosomalCH,
    RecessiveAutosomalHomo,
    XDominant,
    XPseudoDominantFemale,
    XRecessiveFemaleCH,
    XRecessiveFemaleHom,
    XRecessiveMale,
)
from talos.pedigree_parser import PedigreeParser

from test.test_utils import TWO_EXPECTED

TEST_COORDS = Coordinates(chrom='1', pos=1, ref='A', alt='C')
TEST_COORDS2 = Coordinates(chrom='2', pos=2, ref='G', alt='T')
TEST_COORDS_X_1 = Coordinates(chrom='X', pos=1, ref='G', alt='T')
TEST_COORDS_X_2 = Coordinates(chrom='X', pos=2, ref='G', alt='T')

# just to please Ruff
TWO_RESULTS = 2


@pytest.mark.parametrize(
    'moi_string,filters',
    (
        ('Monoallelic', ['DominantAutosomal']),
        ('Mono_And_Biallelic', ['DominantAutosomal', 'RecessiveAutosomalCH', 'RecessiveAutosomalHomo']),
        ('Unknown', ['DominantAutosomal', 'RecessiveAutosomalCH', 'RecessiveAutosomalHomo']),
        ('Biallelic', ['RecessiveAutosomalCH', 'RecessiveAutosomalHomo']),
        ('Hemi_Mono_In_Female', ['XRecessiveMale', 'XDominant']),
        (
            'Hemi_Bi_In_Female',
            ['XRecessiveMale', 'XRecessiveFemaleCH', 'XRecessiveFemaleHom', 'XPseudoDominantFemale'],
        ),
    ),
)
def test_moi_runner(moi_string: str, filters: list[str], pedigree_path):
    """
    check that the right methods are associated with each MOI
    """
    test_runner = MOIRunner(pedigree=PedigreeParser(pedigree_path), target_moi=moi_string)

    # the imported (uninstantiated) objects don't have __class__
    # and the instantiated objects don't have a __name__
    for each_filter in test_runner.filter_list:
        assert each_filter.__class__.__name__ in filters


def test_dominant_autosomal_fails_on_depth(pedigree_path):
    """
    test case for autosomal dominant depth failure
    """

    info_dict = {'gnomad_af': 0.0001, 'gnomad_ac': 0, 'gnomad_hom': 0, 'gene_id': 'TEST1'}
    dom = DominantAutosomal(pedigree=PedigreeParser(pedigree_path))

    # passes with heterozygous
    shallow_variant = SmallVariant(
        info=info_dict,
        het_samples={'male'},
        coordinates=TEST_COORDS,
        depths={'male': 1},
        transcript_consequences=[],
    )
    results = dom.run(principal=shallow_variant)
    assert len(results) == 0


def test_dominant_autosomal_passes(pedigree_path):
    """
    test case for autosomal dominant
    :return:
    """

    info_dict = {
        'gnomad_af': 0.00001,
        'ac': 10,
        'af': 0.0001,
        'gnomad_ac': 0,
        'gnomad_hom': 0,
        'categoryboolean1': True,
        'gene_id': 'TEST1',
    }

    # attributes relating to categorisation
    boolean_categories = ['categoryboolean1']

    dom = DominantAutosomal(pedigree=PedigreeParser(pedigree_path))

    # passes with heterozygous
    passing_variant = SmallVariant(
        info=info_dict,
        het_samples={'male'},
        coordinates=TEST_COORDS,
        boolean_categories=boolean_categories,
        alt_depths={'male': 109},
        depths={'male': 999},
        transcript_consequences=[],
    )
    results = dom.run(principal=passing_variant)
    assert len(results) == 1
    assert results[0].reasons == 'Autosomal Dominant'

    # also passes with homozygous
    passing_variant = SmallVariant(
        info=info_dict,
        hom_samples={'male'},
        coordinates=TEST_COORDS,
        boolean_categories=boolean_categories,
        alt_depths={'male': 199},
        depths={'male': 999},
        transcript_consequences=[],
    )
    results = dom.run(principal=passing_variant)
    assert len(results) == 1
    assert results[0].reasons == 'Autosomal Dominant'

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
        {'gnomad_af': 0.1, 'ac': 10, 'af': 0.0001},
        {'gnomad_hom': 2, 'ac': 10, 'af': 0.0001},
    ],
)
def test_dominant_autosomal_fails(info, pedigree_path):
    """
    test case for autosomal dominant
    :param info: info dict for the variant
    :return:
    """

    dom = DominantAutosomal(pedigree=PedigreeParser(pedigree_path))

    # fails due to high af
    failing_variant = SmallVariant(info=info, het_samples={'male'}, coordinates=TEST_COORDS, transcript_consequences=[])
    assert not dom.run(principal=failing_variant)


def test_recessive_autosomal_hom_passes(pedigree_path):
    """
    check that when the info values are defaults (0)
    we accept a homozygous variant as a Recessive
    """
    passing_variant = SmallVariant(
        hom_samples={'male'},
        coordinates=TEST_COORDS,
        ab_ratios={'male': 1.0},
        alt_depths={'male': 15},
        depths={'male': 15},
        boolean_categories=['categoryboolean1'],
        info={'categoryboolean1': True, 'gene_id': 'TEST1', 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    rec = RecessiveAutosomalHomo(pedigree=PedigreeParser(pedigree_path))
    results = rec.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == 'Autosomal Recessive Homozygous'


def test_recessive_autosomal_hom_fails_alt_reads(pedigree_path):
    """
    check that when the info values are defaults (0)
    we accept a homozygous variant as a Recessive
    """
    passing_variant = SmallVariant(
        hom_samples={'male'},
        coordinates=TEST_COORDS,
        ab_ratios={'male': 1.0},
        alt_depths={'male': 4},
        depths={'male': 15},
        boolean_categories=['categoryboolean1'],
        info={'categoryboolean1': True, 'gene_id': 'TEST1', 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    rec = RecessiveAutosomalHomo(pedigree=PedigreeParser(pedigree_path))
    results = rec.run(passing_variant)
    assert not results


def test_recessive_autosomal_hom_passes_with_ab_flag(pedigree_path):
    """
    check that when the info values are defaults (0)
    we accept a homozygous variant as a Recessive
    """

    passing_variant = SmallVariant(
        hom_samples={'male'},
        coordinates=TEST_COORDS,
        ab_ratios={'male': 0.4},
        alt_depths={'male': 40},
        depths={'male': 40},
        boolean_categories=['categoryboolean1'],
        info={'categoryboolean1': True, 'gene_id': 'TEST1', 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    rec = RecessiveAutosomalHomo(pedigree=PedigreeParser(pedigree_path))
    results = rec.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == 'Autosomal Recessive Homozygous'
    assert passing_variant.get_sample_flags('male') == {'AB Ratio'}


def test_recessive_autosomal_comp_het_male_passes(pedigree_path):
    """
    check that when the info values are defaults (0)
    and the comp-het test is always True
    we accept a heterozygous variant as a Comp-Het
    """

    passing_variant = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS,
        ab_ratios={'male': 0.5},
        alt_depths={'male': 25},
        depths={'male': 50},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True, 'var_link': 'passing1', 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    passing_variant2 = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS2,
        ab_ratios={'male': 0.5},
        alt_depths={'male': 25},
        depths={'male': 50},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True, 'var_link': 'passing2', 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    comp_hets = {'male': {TEST_COORDS.string_format: [passing_variant2]}}
    rec = RecessiveAutosomalCH(pedigree=PedigreeParser(pedigree_path))
    results = rec.run(passing_variant, comp_het=comp_hets)
    assert len(results) == 1
    assert results[0].reasons == 'Autosomal Recessive Comp-Het'


def test_recessive_autosomal_comp_het_male_passes_with_support(pedigree_path):
    """
    check that when one variant is support only, and the other is full, both reach the report
    """

    passing_variant = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS,
        ab_ratios={'male': 0.5},
        alt_depths={'male': 25},
        depths={'male': 50},
        boolean_categories=['categoryboolean1'],
        support_categories={'6'},
        info={'gene_id': 'TEST1', 'categoryboolean1': True, 'var_link': 'passing1', 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    passing_variant2 = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS2,
        ab_ratios={'male': 0.5},
        alt_depths={'male': 25},
        depths={'male': 50},
        boolean_categories=['categoryboolean6'],
        support_categories={'6'},
        info={'gene_id': 'TEST1', 'categoryboolean6': True, 'var_link': 'passing2', 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    comp_hets = {
        'male': {
            TEST_COORDS.string_format: [passing_variant2],
            TEST_COORDS2.string_format: [passing_variant],
        },
    }
    rec = RecessiveAutosomalCH(pedigree=PedigreeParser(pedigree_path))
    results = rec.run(passing_variant, comp_het=comp_hets)
    results.extend(rec.run(passing_variant2, comp_het=comp_hets))
    assert len(results) == TWO_RESULTS
    assert results[0].reasons == 'Autosomal Recessive Comp-Het'


def test_recessive_autosomal_comp_het_male_fails_both_support(pedigree_path):
    """
    check that when one variant is support only, and the other is full, both reach the report
    """

    passing_variant = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS,
        ab_ratios={'male': 0.5},
        alt_depths={'male': 25},
        depths={'male': 50},
        boolean_categories=['categoryboolean6'],
        support_categories={'6'},
        info={'gene_id': 'TEST1', 'categoryboolean6': True, 'var_link': 'passing1', 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    passing_variant2 = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS2,
        ab_ratios={'male': 0.5},
        alt_depths={'male': 25},
        depths={'male': 50},
        boolean_categories=['categoryboolean6'],
        support_categories={'6'},
        info={'gene_id': 'TEST1', 'categoryboolean6': True, 'var_link': 'passing2', 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    comp_hets = {
        'male': {
            TEST_COORDS.string_format: [passing_variant2],
            TEST_COORDS2.string_format: [passing_variant],
        },
    }
    rec = RecessiveAutosomalCH(pedigree=PedigreeParser(pedigree_path))
    results = rec.run(passing_variant, comp_het=comp_hets)
    results.extend(rec.run(passing_variant2, comp_het=comp_hets))
    assert len(results) == 0


def test_recessive_autosomal_comp_het_male_passes_partner_flag(pedigree_path):
    """
    info values are defaults (0) & comp-het test is always True we accept a heterozygous variant as a Comp-Het
    """

    passing_variant = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS,
        ab_ratios={'male': 0.5},
        alt_depths={'male': 25},
        depths={'male': 50},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True, 'var_link': 'passing1', 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    passing_variant2 = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS2,
        ab_ratios={'male': 1.0},
        alt_depths={'male': 25},
        depths={'male': 50},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True, 'var_link': 'passing2', 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    comp_hets = {'male': {TEST_COORDS.string_format: [passing_variant2]}}
    rec = RecessiveAutosomalCH(pedigree=PedigreeParser(pedigree_path))
    results = rec.run(passing_variant, comp_het=comp_hets)
    assert len(results) == 1
    assert results[0].reasons == 'Autosomal Recessive Comp-Het'
    assert results[0].flags == {'AB Ratio'}
    assert results[0].support_vars == {'passing2'}


def test_recessive_autosomal_comp_het_female_passes(pedigree_path):
    """
    info values are defaults (0) & comp-het test is always True we accept a heterozygous variant as a Comp-Het
    """
    passing_variant = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS,
        ab_ratios={'female': 0.5},
        alt_depths={'female': 50},
        depths={'female': 50},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True, 'var_link': 'passing1', 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    passing_variant2 = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS2,
        ab_ratios={'female': 0.5},
        alt_depths={'female': 50},
        depths={'female': 50},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True, 'var_link': 'passing2', 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    comp_hets = {'female': {TEST_COORDS.string_format: [passing_variant2]}}
    rec = RecessiveAutosomalCH(pedigree=PedigreeParser(pedigree_path))
    results = rec.run(passing_variant, comp_het=comp_hets)
    assert len(results) == 1
    assert results[0].reasons == 'Autosomal Recessive Comp-Het'
    assert not results[0].flags
    assert results[0].support_vars == {'passing2'}


def test_recessive_autosomal_comp_het_fails_no_ch_return(pedigree_path):
    """
    check that when the info values are defaults (0) & comp-het test is always False we have no accepted MOI
    """
    failing_variant = SmallVariant(
        info={'gene_id': 'TEST1', 'ac': 10, 'af': 0.0001},
        het_samples={'male'},
        alt_depths={'male': 25},
        depths={'male': 50},
        coordinates=TEST_COORDS,
        transcript_consequences=[],
    )
    rec = RecessiveAutosomalCH(pedigree=PedigreeParser(pedigree_path))
    assert not rec.run(failing_variant)


def test_recessive_autosomal_comp_het_fails_no_paired_call(pedigree_path):
    """
    check that when the info values are defaults (0) & comp-het test is False we have no accepted MOI
    """

    failing_variant = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS,
        ab_ratios={'male': 0.5},
        alt_depths={'male': 25},
        depths={'male': 50},
        info={'gene_id': 'TEST1', 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    failing_variant2 = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS2,
        ab_ratios={'female': 0.5},
        alt_depths={'female': 50},
        depths={'female': 50},
        info={'gene_id': 'TEST1'},
        transcript_consequences=[],
    )

    rec = RecessiveAutosomalCH(pedigree=PedigreeParser(pedigree_path))
    assert not rec.run(
        failing_variant,
        comp_het={'male': {TEST_COORDS2.string_format: [failing_variant2]}},
    )


@pytest.mark.parametrize('info', [{'gnomad_hom': 3, 'gene_id': 'TEST1', 'ac': 10, 'af': 0.0001}])  # threshold is 2
def test_recessive_autosomal_hom_fails(info, pedigree_path):
    """
    check that when the info values are failures we have no confirmed MOI
    """
    failing_variant = SmallVariant(
        info=info,
        het_samples={'male'},
        hom_samples={'male'},
        coordinates=TEST_COORDS,
        transcript_consequences=[],
    )
    rec = RecessiveAutosomalHomo(pedigree=PedigreeParser(pedigree_path))
    assert not rec.run(failing_variant)


def test_x_dominant_female_and_male_het_passes(pedigree_path):
    """
    check that a male and female are accepted as dominant hets
    """
    passing_variant = SmallVariant(
        boolean_categories=['categoryboolean1'],
        info={'gnomad_hemi': 0, 'gene_id': 'TEST1', 'categoryboolean1': True, 'ac': 10, 'af': 0.0001},
        het_samples={'female', 'male'},
        alt_depths={'female': 50, 'male': 50},
        depths={'female': 50, 'male': 50},
        coordinates=TEST_COORDS_X_1,
        transcript_consequences=[],
    )
    x_dom = XDominant(pedigree=PedigreeParser(pedigree_path))
    results = x_dom.run(passing_variant)

    assert len(results) == TWO_EXPECTED
    assert results[0].reasons == 'X_Dominant'


def test_x_dominant_female_hom_passes(pedigree_path):
    """
    check that a female is accepted as a hom
    """
    passing_variant = SmallVariant(
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True, 'ac': 10, 'af': 0.0001},
        hom_samples={'female'},
        alt_depths={'female': 100},
        depths={'female': 100},
        ab_ratios={'female': 0.5},
        coordinates=TEST_COORDS_X_1,
        transcript_consequences=[],
    )
    x_dom = XDominant(pedigree=PedigreeParser(pedigree_path))
    results = x_dom.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == 'X_Dominant'


def test_x_dominant_male_hom_passes(pedigree_path):
    """
    check that a male is accepted as a het
    """
    passing_variant = SmallVariant(
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True, 'ac': 10, 'af': 0.0001},
        hom_samples={'male'},
        alt_depths={'male': 100},
        depths={'male': 100},
        coordinates=TEST_COORDS_X_1,
        transcript_consequences=[],
    )
    x_dom = XDominant(pedigree=PedigreeParser(pedigree_path))
    results = x_dom.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == 'X_Dominant'


@pytest.mark.parametrize(
    'info,wins',
    [
        ({'gnomad_af': 0.00001, 'gene_id': 'TEST1', 'categoryboolean1': True, 'ac': 10, 'af': 0.0001}, 1),
        ({'gnomad_homalt': 2, 'gene_id': 'TEST1', 'categoryboolean1': True, 'ac': 10, 'af': 0.0001}, 1),
        ({'gnomad_ac': 11, 'gene_id': 'TEST1', 'categoryboolean1': False, 'ac': 10, 'af': 0.0001}, 0),
    ],
)
def test_x_dominant_info_fails(info: dict, wins: int, pedigree_path):
    """
    check for info dict exclusions
    """
    passing_variant = SmallVariant(
        info=info,
        hom_samples={'male'},
        coordinates=TEST_COORDS_X_1,
        boolean_categories=['categoryboolean1'],
        transcript_consequences=[],
        alt_depths={'male': 100},
        depths={'male': 100},
    )
    x_dom = XDominant(pedigree=PedigreeParser(pedigree_path))
    assert len(x_dom.run(passing_variant)) == wins


def test_x_dominant_female_inactivation_passes(pedigree_path):
    """
    check that a het female, but not a het male, is accepted as dominant
    current test implementation doesn't include family consideration
    """
    info_dict = {
        'gnomad_af': 0.0,
        'gnomad_ac': 0,
        'gnomad_homalt': 0,
        'gene_id': 'TEST1',
        'categoryboolean1': True,
        'ac': 10,
        'af': 0.0001,
    }
    passing_variant = SmallVariant(
        alt_depths={'female': 100, 'male': 100},
        depths={'female': 100, 'male': 100},
        info=info_dict,
        het_samples={'female', 'male'},
        hom_samples={'male'},
        boolean_categories=['categoryboolean1'],
        coordinates=TEST_COORDS_X_1,
        transcript_consequences=[],
    )
    x_dom = XPseudoDominantFemale(pedigree=PedigreeParser(pedigree_path))
    results = x_dom.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == 'X_PseudoDominant'
    assert 'Affected female with heterozygous variant in XLR gene' in results[0].flags


def test_x_recessive_male_hom_passes(pedigree_path):
    passing_variant = SmallVariant(
        hom_samples={'female', 'male'},
        coordinates=TEST_COORDS_X_1,
        ab_ratios={'female': 1.0, 'male': 1.0},
        alt_depths={'female': 100, 'male': 100},
        depths={'female': 100, 'male': 100},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True, 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    x_rec = XRecessiveMale(pedigree=PedigreeParser(pedigree_path))
    results = x_rec.run(passing_variant, comp_het={})
    assert len(results) == 1
    assert results[0].reasons == 'X_Male'


def test_x_recessive_female_hom_passes(pedigree_path):
    """
    :return:
    """

    passing_variant = SmallVariant(
        hom_samples={'female', 'male'},
        coordinates=TEST_COORDS_X_1,
        ab_ratios={'female': 1.0, 'male': 1.0},
        alt_depths={'female': 100, 'male': 100},
        depths={'female': 100, 'male': 100},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True, 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    x_rec = XRecessiveFemaleHom(pedigree=PedigreeParser(pedigree_path))
    results = x_rec.run(passing_variant, comp_het={})
    assert len(results) == 1
    assert results[0].reasons == 'X_Recessive HOM Female'


def test_x_recessive_male_het_passes(pedigree_path):
    passing_variant = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS_X_1,
        ab_ratios={'male': 0.5},
        alt_depths={'male': 25},
        depths={'male': 50},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True, 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    x_rec = XRecessiveMale(pedigree=PedigreeParser(pedigree_path))
    results = x_rec.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == 'X_Male'


def test_x_recessive_female_het_passes(pedigree_path):
    passing_variant = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_1,
        ab_ratios={'female': 0.5},
        alt_depths={'female': 25},
        depths={'female': 50},
        sample_categories=['categorysample4'],
        info={'gene_id': 'TEST1', 'categorysample4': ['female'], 'var_link': 'passing1', 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    passing_variant_2 = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_2,
        ab_ratios={'female': 0.5},
        alt_depths={'female': 50},
        depths={'female': 50},
        sample_categories=['categorysample4'],
        info={'gene_id': 'TEST1', 'categorysample4': ['female'], 'var_link': 'passing2', 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    comp_hets = {'female': {'X-1-G-T': [passing_variant_2]}}
    x_rec = XRecessiveFemaleCH(pedigree=PedigreeParser(pedigree_path))
    results = x_rec.run(passing_variant, comp_het=comp_hets)
    assert len(results) == 1
    assert results[0].reasons == 'X_RecessiveFemaleCompHet'
    assert results[0].support_vars == {'passing2'}


def test_x_recessive_female_het_passes_one_support(pedigree_path):
    passing_variant = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_1,
        ab_ratios={'female': 0.5},
        alt_depths={'female': 25},
        depths={'female': 50},
        sample_categories=['categorysample4'],
        support_categories={'4'},
        info={'gene_id': 'TEST1', 'categorysample4': ['female'], 'var_link': 'passing1', 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    passing_variant_2 = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_2,
        ab_ratios={'female': 0.5},
        alt_depths={'female': 50},
        depths={'female': 50},
        sample_categories=['categorysample4'],
        support_categories=set(),
        info={'gene_id': 'TEST1', 'categorysample4': ['female'], 'var_link': 'passing2', 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    comp_hets = {
        'female': {
            TEST_COORDS_X_1.string_format: [passing_variant_2],
            TEST_COORDS_X_2.string_format: [passing_variant],
        },
    }
    x_rec = XRecessiveFemaleCH(pedigree=PedigreeParser(pedigree_path))
    results = x_rec.run(passing_variant, comp_het=comp_hets)
    results.extend(x_rec.run(passing_variant_2, comp_het=comp_hets))
    assert len(results) == TWO_RESULTS
    assert results[0].reasons == 'X_RecessiveFemaleCompHet'
    assert results[0].support_vars == {'passing2'}


def test_x_recessive_female_het_fails_both_support(pedigree_path):
    passing_variant = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_1,
        ab_ratios={'female': 0.5},
        alt_depths={'female': 25},
        depths={'female': 50},
        sample_categories=['categorysample4'],
        support_categories={'4'},
        info={'gene_id': 'TEST1', 'categorysample4': ['female'], 'var_link': 'passing1', 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    passing_variant_2 = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_2,
        ab_ratios={'female': 0.5},
        alt_depths={'female': 50},
        depths={'female': 50},
        sample_categories=['categorysample4'],
        support_categories={'4'},
        info={'gene_id': 'TEST1', 'categorysample4': ['female'], 'var_link': 'passing2', 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    comp_hets = {
        'female': {
            TEST_COORDS_X_1.string_format: [passing_variant_2],
            TEST_COORDS_X_2.string_format: [passing_variant],
        },
    }
    x_rec = XRecessiveFemaleCH(pedigree=PedigreeParser(pedigree_path))
    results = x_rec.run(passing_variant, comp_het=comp_hets)
    results.extend(x_rec.run(passing_variant_2, comp_het=comp_hets))
    assert len(results) == 0


def test_het_de_novo_passes(pedigree_path):
    passing_variant = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_1,
        sample_categories=['categorysample4'],
        ab_ratios={'female': 0.5},
        alt_depths={'female': 40},
        depths={'female': 99},
        info={'gene_id': 'TEST1', 'categorysample4': ['female'], 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    dom_a = DominantAutosomal(pedigree=PedigreeParser(pedigree_path))
    results = dom_a.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == 'Autosomal Dominant'
    assert not results[0].flags


def test_het_de_novo_het_passes_flagged(pedigree_path):
    passing_variant = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_1,
        sample_categories=['categorysample4'],
        ab_ratios={'female': 0.5},
        alt_depths={'female': 40},
        depths={'female': 99},
        info={'gene_id': 'TEST1', 'categorysample4': ['female'], 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    dom_a = DominantAutosomal(pedigree=PedigreeParser(pedigree_path))
    results = dom_a.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == 'Autosomal Dominant'


def test_x_recessive_female_het_fails(pedigree_path):
    passing_variant = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_1,
        ab_ratios={'female': 0.5},
        alt_depths={'female': 25},
        depths={'female': 50},
        sample_categories=['categorysample4'],
        info={'gene_id': 'TEST1', 'categorysample4': ['male'], 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    passing_variant_2 = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS_X_2,
        ab_ratios={'male': 0.5},
        alt_depths={'male': 25},
        depths={'male': 50},
        sample_categories=['categorysample4'],
        info={'gene_id': 'TEST1', 'categorysample4': ['male'], 'ac': 10, 'af': 0.0001},
        transcript_consequences=[],
    )
    comp_hets = {'female': {'x-2-A-C': [passing_variant_2]}}
    x_rec = XRecessiveFemaleCH(pedigree=PedigreeParser(pedigree_path))
    assert not x_rec.run(passing_variant, comp_het=comp_hets)


def test_x_recessive_female_het_no_pair_fails(pedigree_path):
    """ """

    passing_variant = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_1,
        ab_ratios={'female': 0.5},
        alt_depths={'female': 25},
        depths={'female': 50},
        info={
            'gene_id': 'TEST1',
            'categorysample1': True,
            'boolean_categories': 'categorysample1',
            'ac': 10,
            'af': 0.0001,
        },
        transcript_consequences=[],
    )
    assert not XRecessiveFemaleCH(pedigree=PedigreeParser(pedigree_path)).run(passing_variant)


def test_check_familial_inheritance_simple(pedigree_path):
    """
    test the check_familial_inheritance method
    trio male, mother_1, father_1; only 'male' is affected
    """

    base_moi = BaseMoi(pedigree=PedigreeParser(pedigree_path), applied_moi='applied')
    assert base_moi.single_variant_explains_disease_in_family(sample_id='male', called_variants={'male'})


def test_check_familial_inheritance_mother_fail(pedigree_path):
    """
    test the check_familial_inheritance method
    """

    base_moi = BaseMoi(pedigree=PedigreeParser(pedigree_path), applied_moi='applied')
    assert not base_moi.single_variant_explains_disease_in_family(
        sample_id='male',
        called_variants={'male', 'mother_1'},
    )


def test_check_familial_inheritance_mother_passes(pedigree_path):
    """
    test the check_familial_inheritance method
    mother in variant calls, but partial penetrance
    """

    base_moi = BaseMoi(pedigree=PedigreeParser(pedigree_path), applied_moi='applied')

    assert base_moi.single_variant_explains_disease_in_family(
        sample_id='male',
        called_variants={'male', 'mother_1'},
        partial_pen=True,
    )


def test_check_familial_inheritance_father_fail(pedigree_path):
    """
    test the check_familial_inheritance method
    """

    base_moi = BaseMoi(pedigree=PedigreeParser(pedigree_path), applied_moi='applied')
    assert not base_moi.single_variant_explains_disease_in_family(
        sample_id='male',
        called_variants={'male', 'father_1'},
    )


def test_check_familial_inheritance_father_passes(pedigree_path):
    """
    test the check_familial_inheritance method
    father in variant calls, but partial penetrance
    """

    base_moi = BaseMoi(pedigree=PedigreeParser(pedigree_path), applied_moi='applied')

    result = base_moi.single_variant_explains_disease_in_family(
        sample_id='male',
        called_variants={'male', 'father_1'},
        partial_pen=True,
    )
    assert result


def test_check_familial_inheritance_top_down(pedigree_path):
    """
    test the check_familial_inheritance method
    father in variant calls, but partial penetrance
    """

    base_moi = BaseMoi(pedigree=PedigreeParser(pedigree_path), applied_moi='applied')
    assert base_moi.single_variant_explains_disease_in_family(
        sample_id='father_1',
        called_variants={'male', 'father_1'},
        partial_pen=True,
    )


def test_check_familial_inheritance_no_calls(pedigree_path):
    """
    test the check_familial_inheritance method where there are no calls
    we lazily pass this - the assumption is that we're only assessing samples with variants
    this method checks parents as a trio, not the sample in question
    """

    base_moi = BaseMoi(pedigree=PedigreeParser(pedigree_path), applied_moi='applied')
    assert base_moi.single_variant_explains_disease_in_family(sample_id='male', called_variants=set(), partial_pen=True)


def test_genotype_calls(pedigree_path):
    """
    test the manual genotype assignments
    """
    base_moi = DominantAutosomal(pedigree=PedigreeParser(pedigree_path), applied_moi='applied')

    info_dict = {'gnomad_af': 0.0001, 'gnomad_ac': 0, 'gnomad_hom': 0, 'gene_id': 'TEST1'}
    variant = SmallVariant(
        info=info_dict,
        het_samples={'male'},
        hom_samples={'female'},
        coordinates=TEST_COORDS,
        transcript_consequences=[],
    )
    assert base_moi.get_family_genotypes(variant, 'male') == {'father_1': 'WT', 'male': 'Het', 'mother_1': 'WT'}
    assert base_moi.get_family_genotypes(variant, 'female') == {'father_2': 'WT', 'female': 'Hom', 'mother_2': 'WT'}
    x_variant = SmallVariant(
        info=info_dict,
        het_samples={'male', 'female'},
        coordinates=TEST_COORDS_X_1,
        transcript_consequences=[],
    )
    assert base_moi.get_family_genotypes(x_variant, 'male') == {'father_1': 'WT', 'male': 'Hemi', 'mother_1': 'WT'}
    assert base_moi.get_family_genotypes(x_variant, 'female') == {'father_2': 'WT', 'female': 'Het', 'mother_2': 'WT'}

    x_variant_2 = SmallVariant(
        info=info_dict,
        hom_samples={'male', 'female'},
        coordinates=TEST_COORDS_X_1,
        transcript_consequences=[],
    )
    assert base_moi.get_family_genotypes(x_variant_2, 'male') == {'father_1': 'WT', 'male': 'Hemi', 'mother_1': 'WT'}
    assert base_moi.get_family_genotypes(x_variant_2, 'female') == {'father_2': 'WT', 'female': 'Hom', 'mother_2': 'WT'}
    variant_missing = SmallVariant(info=info_dict, coordinates=TEST_COORDS, transcript_consequences=[])
    assert base_moi.get_family_genotypes(variant_missing, 'male') == {'father_1': 'WT', 'male': 'WT', 'mother_1': 'WT'}


def test_too_common_in_population_passes_cat1():
    """this should return false, i.e. would not be removed"""
    info = {'categoryboolean1': True, 'ac': 10, 'af': 0.0001}
    assert not ClinVarFilter().too_common(
        variant=SmallVariant(
            info=info,
            coordinates=TEST_COORDS,
            transcript_consequences=[],
        ),
    )


def test_too_common_in_population_passes_below_thresholds():
    """should return false, i.e. attributes are below thresholds"""
    info = {'gnomad_af': 0.0, 'ac': 10, 'af': 0.0001}
    assert not ClinVarFilter().too_common(
        variant=SmallVariant(
            info=info,
            coordinates=TEST_COORDS,
            transcript_consequences=[],
        ),
    )


def test_too_common_in_population_true():
    """should return true, i.e. attributes above threshold"""
    info = {'gnomad_af': 0.1}
    assert ClinVarFilter().too_common(
        variant=SmallVariant(
            info=info,
            coordinates=TEST_COORDS,
            transcript_consequences=[],
        ),
    )


def test_too_common_in_population_true_multi():
    """should return true, i.e. attributes above thresholds, even if others are below"""
    info = {'gnomad_af': 0.0, 'ac': 11, 'af': 0.1}
    assert ClinVarFilter().too_common(
        variant=SmallVariant(
            info=info,
            coordinates=TEST_COORDS,
            transcript_consequences=[],
        ),
    )
