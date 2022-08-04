"""
test class for the utils collection
"""
from dataclasses import dataclass
from typing import List
import pytest
from cyvcf2 import VCFReader
from reanalysis.utils import (
    AbstractVariant,
    find_comp_hets,
    gather_gene_dict_from_contig,
    get_non_ref_samples,
    get_simple_moi,
    read_json_from_path,
    identify_file_type,
    FileTypes,
)


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
    with pytest.raises(Exception):
        identify_file_type('i/am/a/mystery.file.type')


def test_read_json(conf_json_path):
    """

    :return:
    """
    parsed_json = read_json_from_path(conf_json_path)
    assert 'moi_tests' in parsed_json.keys()

    with pytest.raises(FileNotFoundError):
        read_json_from_path('not_a_real_path')


@pytest.mark.parametrize(
    'string,expected',
    [
        ('blag', 'Mono_And_Biallelic'),
        (None, 'Mono_And_Biallelic'),
        ('Unknown', 'Mono_And_Biallelic'),
        ('BIALLELIC_ANY', 'Biallelic'),
        ('BOTH_something,something', 'Mono_And_Biallelic'),
        ('MONO,something', 'Monoallelic'),
        ('X-LINKED', 'Hemi_Mono_In_Female'),
        ('X-LINKED biallelic', 'Hemi_Bi_In_Female'),
    ],
)
def test_get_simple_moi(string: str, expected: str):
    """

    :param string:
    :param expected:
    :return:
    """
    assert get_simple_moi(string) == expected


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


def test_gene_dict(two_trio_variants_vcf, conf_json_path):
    """
    gene = ENSG00000075043
    :param two_trio_variants_vcf:
    :param conf_json_path:
    :return:
    """
    conf = read_json_from_path(conf_json_path)
    reader = VCFReader(two_trio_variants_vcf)
    contig = 'chr20'
    panel_data = {'ENSG00000075043': {'new': True}}
    var_dict = gather_gene_dict_from_contig(
        contig=contig, config=conf, variant_source=reader, panelapp_data=panel_data
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


def test_phased_dict(phased_vcf_path, conf_json_path):
    """
    gene = ENSG00000075043
    :param phased_vcf_path:
    :param conf_json_path:
    :return:
    """
    conf = read_json_from_path(conf_json_path)
    reader = VCFReader(phased_vcf_path)
    panel_data = {'ENSG00000075043': {'new': True}}
    var_dict = gather_gene_dict_from_contig(
        contig='chr20', config=conf, variant_source=reader, panelapp_data=panel_data
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
