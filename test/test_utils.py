"""
test class for the utils collection
"""
from dataclasses import dataclass
from typing import List
import os
import pytest

from reanalysis.utils import (
    AbstractVariant,
    get_non_ref_samples,
    get_simple_moi,
    read_json_from_path,
)


PWD = os.path.dirname(__file__)
INPUT = os.path.join(PWD, 'input')
JSON_STUB = os.path.join(INPUT, 'json_example.json')


def test_read_json():
    """

    :return:
    """
    assert read_json_from_path(JSON_STUB) == {'key': 'value'}
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
    assert trio_abs_variant.sample_de_novo('PROBAND')
    assert not trio_abs_variant.sample_de_novo('FATHER')


def test_av_categories(trio_abs_variant: AbstractVariant):
    """

    :param trio_abs_variant:
    :return:
    """
    assert trio_abs_variant.category_4  # non-empty list
    assert trio_abs_variant.category_3
    assert not trio_abs_variant.category_2
    assert not trio_abs_variant.category_1
    assert trio_abs_variant.category_1_2_3
    assert trio_abs_variant.category_non_support
    assert trio_abs_variant.is_classified
