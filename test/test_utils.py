"""
test class for the utils collection
"""

from copy import deepcopy

import pytest
from cyvcf2 import VCFReader

from talos.models import (
    Coordinates,
    FileTypes,
    PanelApp,
    PhenotypeMatchedPanels,
    ReportVariant,
    SmallVariant,
)
from talos.utils import get_simple_moi
from talos.utils import (
    find_comp_hets,
    gather_gene_dict_from_contig,
    get_new_gene_map,
    get_non_ref_samples,
    identify_file_type,
    make_flexible_pedigree,
)

ZERO_EXPECTED = 0
ONE_EXPECTED = 1
TWO_EXPECTED = 2
THREE_EXPECTED = 3
FOUR_EXPECTED = 4
FIVE_EXPECTED = 5


def test_coord_sorting():
    """
    check that coord sorting methods work
    """
    coord_1 = Coordinates(chrom='4', pos=20, ref='A', alt='C')
    coord_1b = Coordinates(chrom='4', pos=21, ref='A', alt='C')
    coord_1c = Coordinates(chrom='4', pos=21, ref='A', alt='C')
    coord_2 = Coordinates(chrom='5', pos=20, ref='A', alt='C')
    assert coord_1 < coord_2
    assert coord_1 < coord_1b
    assert not coord_1b < coord_1c


def test_abs_var_sorting(two_trio_abs_variants: list[SmallVariant]):
    """
    test sorting and equivalence at the AbsVar level
    """

    var1, var2 = two_trio_abs_variants
    assert var1 < var2
    assert sorted([var2, var1]) == [var1, var2]
    # not sure if I should be able to just override the chrom...
    var1.coordinates.chrom = 'HLA1234'
    assert var1 > var2


def test_reported_variant_ordering(trio_abs_variant: SmallVariant):
    """
    test that equivalence between Report objects works as exp.
    """
    report_1 = ReportVariant(
        sample='1',
        family='1',
        gene='2',
        var_data=deepcopy(trio_abs_variant),
        reasons={'test'},
        genotypes={},
    )
    report_2 = ReportVariant(
        sample='1',
        family='1',
        gene='2',
        var_data=deepcopy(trio_abs_variant),
        reasons={'test'},
        genotypes={},
    )
    assert report_1 == report_2
    # alter sample ID, expected mismatch
    report_1.sample = '2'
    assert report_1 != report_2
    report_2.sample = '2'
    report_1.var_data.coordinates.chrom = '1'
    report_2.var_data.coordinates.chrom = '11'
    assert report_1 < report_2


def test_file_types():
    """
    check 'em
    """
    assert identify_file_type('this/is/my/matrixtable.mt') == FileTypes.MATRIX_TABLE
    assert identify_file_type('this/is/my/hailtable.ht') == FileTypes.HAIL_TABLE
    assert identify_file_type('this/is/a/varfile.vcf') == FileTypes.VCF
    assert identify_file_type('this/is/a/varfile.vcf.gz') == FileTypes.VCF_GZ
    assert identify_file_type('this/is/a/varfile.vcf.bgz') == FileTypes.VCF_BGZ


def test_file_types_assert_error():
    """
    check 'em
    """
    with pytest.raises(ValueError):
        identify_file_type('no/extensions')


def test_file_types_exception():
    """
    check 'em
    :return:
    """
    with pytest.raises(TypeError):
        identify_file_type('i/am/a/mystery.file.type')


@pytest.mark.parametrize(
    'string,expected,chrom',
    [
        ('blag', 'Biallelic', '1'),
        ('blag', 'Hemi_Bi_In_Female', 'X'),
        ('biallelic ANY', 'Biallelic', '1'),
        ('both something,something', 'Mono_And_Biallelic', '1'),
        ('monoallelic, something', 'Monoallelic', '1'),
        ('x-linked', 'Hemi_Mono_In_Female', 'X'),
        ('x-linked biallelic', 'Hemi_Bi_In_Female', 'X'),
    ],
)
def test_get_simple_moi(string: str, expected: str, chrom: str):
    """
    Tests the string parsing down to simple representation
    """
    assert get_simple_moi({string}, chrom) == {expected}


def test_get_non_ref_samples(cyvcf_example_variant):
    """
    this simple test can be done without the use of a cyvcf2 object
    :return:
    """

    samples = ['male', 'father', 'mother']
    het, hom = get_non_ref_samples(variant=cyvcf_example_variant, samples=samples)
    assert het == {'male'}
    assert not hom


def test_av_categories(trio_abs_variant: SmallVariant):
    """
    Cat. 3, and Cat. 4 for PROBAND only:
    """
    assert trio_abs_variant.info.get('categoryboolean3')
    assert not trio_abs_variant.info.get('categoryboolean1')
    assert not trio_abs_variant.info.get('categoryboolean2')
    assert trio_abs_variant.sample_category_check('male')

    for sample_cat in trio_abs_variant.sample_categories:
        assert 'father_1' not in trio_abs_variant.info[sample_cat]


def test_av_categories_support(trio_abs_variant: SmallVariant):
    """
    Cat. 3, and Cat. 4 for PROBAND only
    """
    assert trio_abs_variant.info.get('categoryboolean3')
    assert trio_abs_variant.sample_category_check('male')

    # now make the categories support-only
    trio_abs_variant.support_categories.update({'3', '4'})
    assert trio_abs_variant.sample_category_check('male')
    assert not trio_abs_variant.sample_category_check('male', allow_support=False)


def test_av_phase(trio_abs_variant: SmallVariant):
    """
    nothing here yet
    """
    assert trio_abs_variant.phased == {}


def test_gene_dict(two_trio_variants_vcf):
    """
    gene = ENSG00000075043
    """
    reader = VCFReader(two_trio_variants_vcf)
    contig = 'chr20'
    var_dict = gather_gene_dict_from_contig(contig=contig, variant_source=reader)
    assert len(var_dict) == 1
    assert 'ENSG00000075043' in var_dict
    assert len(var_dict['ENSG00000075043']) == TWO_EXPECTED


def test_comp_hets(two_trio_abs_variants: list[SmallVariant], pedigree_path):
    """
    {
        'male': {
            '20-63406931-C-CGG': [Variant()],
            '20-63406991-C-CGG': [Variant()]
        }
    }
    :param two_trio_abs_variants:
    :return:
    """
    ch_dict = find_comp_hets(two_trio_abs_variants, pedigree=make_flexible_pedigree(pedigree_path))
    assert 'male' in ch_dict
    results = ch_dict.get('male')
    assert isinstance(results, dict)
    assert len(results) == TWO_EXPECTED
    key_1, key_2 = list(results.keys())
    assert results[key_1][0].coordinates.string_format == key_2
    assert results[key_2][0].coordinates.string_format == key_1


def test_phased_dict(phased_vcf_path):
    """
    gene = ENSG00000075043
    """
    reader = VCFReader(phased_vcf_path)
    var_dict = gather_gene_dict_from_contig(contig='chr20', variant_source=reader)
    assert len(var_dict) == ONE_EXPECTED
    assert 'ENSG00000075043' in var_dict
    assert len(var_dict['ENSG00000075043']) == TWO_EXPECTED
    var_pair = var_dict['ENSG00000075043']
    for variant in var_pair:
        assert 'mother_1' in variant.phased
        assert variant.phased['mother_1'] == {420: '0|1'}


def test_phased_comp_hets(phased_variants: list[SmallVariant], pedigree_path: str):
    """
    phased variants shouldn't form a comp-het
    'mother_1' is het for both variants, but phase-set is same for both
    :param phased_variants:
    :return:
    """
    ch_dict = find_comp_hets(phased_variants, pedigree=make_flexible_pedigree(pedigree_path))
    assert len(ch_dict) == ZERO_EXPECTED


def test_new_gene_map_null():
    """
    with no specific pheno data, new at all is new for all
    """

    panel_data = PanelApp(genes={'ENSG1': {'new': {1, 2}, 'symbol': 'ensg1'}})
    result = get_new_gene_map(panel_data)
    assert result == {'ENSG1': {'all'}}


def test_new_gene_map_core():
    """
    for a core panel this should also be new
    even if the core panel isn't assigned to individuals
    """

    panel_data = PanelApp(genes={'ENSG1': {'new': {137}, 'symbol': 'ensg1'}})
    personal_panels = PhenotypeMatchedPanels()
    result = get_new_gene_map(panel_data, personal_panels)
    assert result == {'ENSG1': {'all'}}


def test_new_gene_map_cohort_level():
    """
    check that new for the cohort-matched panel is new for all
    """

    panel_data = PanelApp(genes={'ENSG1': {'new': {99}, 'symbol': 'ensg1'}})
    personal_panels = PhenotypeMatchedPanels()
    result = get_new_gene_map(panel_data, personal_panels)
    assert result == {'ENSG1': {'all'}}


def test_new_gene_map_mix_n_match():
    """
    now test the pheno-matched new
    """
    panel_data = PanelApp(genes={'ENSG1': {'new': {1}, 'symbol': 'ensg1'}})
    personal_panels = PhenotypeMatchedPanels(samples={'sam': {'panels': {1, 2}}})
    result = get_new_gene_map(panel_data, personal_panels)
    assert result == {'ENSG1': {'sam'}}


def test_new_gene_map_fail_handled():
    """
    What if we find a panel that wasn't assigned to anyone
    """
    panel_data = PanelApp(genes={'ENSG1': {'new': {2}, 'symbol': 'ensg1'}})
    personal_panels = PhenotypeMatchedPanels(samples={'sam': {'panels': {1}}})
    with pytest.raises(AssertionError):
        get_new_gene_map(panel_data, personal_panels)


def test_new_gene_map_complex():
    """
    ENSG2 is new for everyone
    """
    panel_data = PanelApp(
        genes={
            'ENSG1': {'new': {1}, 'symbol': 'ensg1'},
            'ENSG2': {'new': {137}, 'symbol': 'ensg2'},
            'ENSG3': {'new': {4}, 'symbol': 'ensg3'},
            'ENSG4': {'new': {2}, 'symbol': 'ensg4'},
        },
    )
    personal_panels = PhenotypeMatchedPanels(samples={'sam': {'panels': {1, 2}}, 'sam2': {'panels': {4, 2}}})
    result = get_new_gene_map(panel_data, personal_panels)
    assert result == {
        'ENSG1': {'sam'},
        'ENSG2': {'all'},
        'ENSG3': {'sam2'},
        'ENSG4': {'sam', 'sam2'},
    }
