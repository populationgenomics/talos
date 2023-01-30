"""
test class for the utils collection
"""
from copy import deepcopy
from dataclasses import dataclass
from typing import List
import pytest
from cyvcf2 import VCFReader
from reanalysis.utils import (
    AbstractVariant,
    Coordinates,
    find_comp_hets,
    gather_gene_dict_from_contig,
    get_non_ref_samples,
    get_simple_moi,
    identify_file_type,
    FileTypes,
    ReportedVariant,
)


def test_coord_sorting():
    """
    check that coord sorting methods work
    """
    coord_1 = Coordinates('4', 20, 'A', 'C')
    coord_1b = Coordinates('4', 21, 'A', 'C')
    coord_1c = Coordinates('4', 21, 'A', 'C')
    coord_2 = Coordinates('5', 20, 'A', 'C')
    assert coord_1 < coord_2
    assert coord_1 < coord_1b
    assert not coord_1b < coord_1c  # pylint: disable=unneeded-not


def test_abs_var_sorting(two_trio_abs_variants: list[AbstractVariant]):
    """
    test sorting and equivalence at the AbsVar level
    """

    var1, var2 = two_trio_abs_variants
    assert var1 < var2
    assert var1 == var1  # pylint: disable=comparison-with-itself
    assert sorted([var2, var1]) == [var1, var2]
    # not sure if I should be able to just override the chrom...
    var1.coords.chrom = 'HLA1234'
    assert var1 > var2


def test_reported_variant_ordering(trio_abs_variant):
    """
    test that equivalence between Report objects works as exp.
    """
    report_1 = ReportedVariant(
        sample='1',
        family='1',
        gene='2',
        var_data=deepcopy(trio_abs_variant),
        reasons={'test'},
        supported=False,
        genotypes={},
    )
    report_2 = ReportedVariant(
        sample='1',
        family='1',
        gene='2',
        var_data=deepcopy(trio_abs_variant),
        reasons={'test'},
        supported=False,
        genotypes={},
    )
    assert report_1 == report_2
    report_1.support_vars = ['support']
    assert report_1 != report_2
    # set to same and re-check
    report_2.support_vars = ['support']
    assert report_1 == report_2
    # alter sample ID, expected mismatch
    report_1.sample = '2'
    assert report_1 != report_2
    report_2.sample = '2'
    report_1.var_data.coords.chrom = '1'
    report_2.var_data.coords.chrom = '11'
    assert report_1 < report_2


def test_file_types():
    """
    check 'em
    :return:
    """
    assert identify_file_type('this/is/my/matrixtable.mt') == FileTypes.MATRIX_TABLE
    assert identify_file_type('this/is/my/hailtable.ht') == FileTypes.HAIL_TABLE
    assert identify_file_type('this/is/a/varfile.vcf') == FileTypes.VCF
    assert identify_file_type('this/is/a/varfile.vcf.gz') == FileTypes.VCF_GZ
    assert identify_file_type('this/is/a/varfile.vcf.bgz') == FileTypes.VCF_BGZ


def test_file_types_assert_error():
    """
    check 'em
    :return:
    """
    with pytest.raises(AssertionError):
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
        (None, 'Biallelic', '1'),
        ('monoallelic, something', 'Monoallelic', '1'),
        ('x-linked', 'Hemi_Bi_In_Female', 'X'),
        (None, 'Hemi_Bi_In_Female', 'X'),
        ('x-linked biallelic', 'Hemi_Bi_In_Female', 'X'),
    ],
)
def test_get_simple_moi(string: str, expected: str, chrom: str):
    """
    Tests the string parsing down to simple representation
    :param string:
    :param expected:
    :param chrom:
    """
    assert get_simple_moi(string, chrom) == expected


def test_get_non_ref_samples():
    """
    this simple test can be done without the use of a cyvcf2 object
    :return:
    """

    @dataclass
    class SuperSimple:
        """test_fixture"""

        gt_types: List[int]

    samples = ['a', 'b', 'c', 'd', 'e']
    variant = SuperSimple([0, 1, 2, 3, 1])
    het, hom = get_non_ref_samples(variant=variant, samples=samples)
    assert het == {'b', 'e'}
    assert hom == {'d'}


def test_av_de_novo(trio_abs_variant: AbstractVariant):
    """

    :param trio_abs_variant:
    :return:
    """
    assert trio_abs_variant.sample_de_novo('male')
    assert not trio_abs_variant.sample_de_novo('father_1')


def test_av_categories(trio_abs_variant: AbstractVariant):
    """

    :param trio_abs_variant:
    :return:
    """
    assert trio_abs_variant.is_classified
    assert trio_abs_variant.category_non_support
    assert trio_abs_variant.has_sample_categories
    assert trio_abs_variant.sample_de_novo('male')
    assert trio_abs_variant.info.get('categoryboolean3')
    assert not trio_abs_variant.info.get('categoryboolean1')
    assert not trio_abs_variant.info.get('categoryboolean2')


def test_av_phase(trio_abs_variant: AbstractVariant):
    """
    nothing here yet
    :param trio_abs_variant:
    :return:
    """
    assert trio_abs_variant.phased == {}


def test_gene_dict(two_trio_variants_vcf):
    """
    gene = ENSG00000075043
    :param two_trio_variants_vcf:
    :return:
    """
    reader = VCFReader(two_trio_variants_vcf)
    contig = 'chr20'
    panel_data = {'ENSG00000075043': {'new': True}}
    var_dict = gather_gene_dict_from_contig(
        contig=contig, variant_source=reader, panelapp_data=panel_data
    )
    assert len(var_dict) == 1
    assert 'ENSG00000075043' in var_dict
    assert len(var_dict['ENSG00000075043']) == 2


def test_comp_hets(two_trio_abs_variants: list[AbstractVariant], peddy_ped):
    """
    {
        'male': {
            '20-63406931-C-CGG': [AbstractVariant()],
            '20-63406991-C-CGG': [AbstractVariant()]
        }
    }
    :param two_trio_abs_variants:
    :return:
    """
    ch_dict = find_comp_hets(two_trio_abs_variants, pedigree=peddy_ped)
    assert 'male' in ch_dict
    results = ch_dict.get('male')
    assert len(results) == 2
    key_1, key_2 = list(results.keys())
    assert results[key_1][0].coords.string_format == key_2
    assert results[key_2][0].coords.string_format == key_1


def test_phased_dict(phased_vcf_path):
    """
    gene = ENSG00000075043
    :param phased_vcf_path:
    :return:
    """
    reader = VCFReader(phased_vcf_path)
    panel_data = {'ENSG00000075043': {'new': True}}
    var_dict = gather_gene_dict_from_contig(
        contig='chr20', variant_source=reader, panelapp_data=panel_data
    )
    assert len(var_dict) == 1
    assert 'ENSG00000075043' in var_dict
    assert len(var_dict['ENSG00000075043']) == 2
    var_pair = var_dict['ENSG00000075043']
    for variant in var_pair:
        assert 'mother_1' in variant.phased
        assert variant.phased['mother_1'] == {420: '0|1'}


def test_phased_comp_hets(phased_variants: list[AbstractVariant], peddy_ped):
    """
    phased variants shouldn't form a comp-het
    'mother_1' is het for both variants, but phase-set is same for both
    :param phased_variants:
    :return:
    """
    ch_dict = find_comp_hets(phased_variants, pedigree=peddy_ped)
    assert len(ch_dict) == 0
