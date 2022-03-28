"""
test class for the utils collection
"""
from dataclasses import dataclass
from typing import List
import os
import pytest

# import cyvcf2

from reanalysis.utils import (
    # AnalysisVariant,
    get_non_ref_samples,
    get_simple_moi,
    read_json_dict_from_path,
    parse_ped_simple,
    PedPerson,
)


PWD = os.path.dirname(__file__)
INPUT = os.path.join(PWD, 'input')
JSON_STUB = os.path.join(INPUT, 'json_example.json')
COMP_HET = os.path.join(INPUT, 'comp_het.vcf.bgz')
SINGLE_VAR = os.path.join(INPUT, 'single_hail.vcf.bgz')
PANELAPP_LATEST = os.path.join(INPUT, 'panelapp_current_137.json')


def test_read_json():
    """

    :return:
    """
    assert read_json_dict_from_path(JSON_STUB) == {'key': 'value'}
    with pytest.raises(FileNotFoundError):
        read_json_dict_from_path('not_a_real_path')


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


def test_pedperson():
    """
    ummm... no methods here, not much to test
    :return:
    """
    test_person = PedPerson('sam', True, False)
    assert test_person.sample == 'sam'
    assert test_person.male
    assert not test_person.affected


def test_parse_ped_simple(tmp_path):
    """
    write then read
    :param tmp_path:
    :return:
    """
    tsv = tmp_path / 'test.tsv'
    with open(tsv, 'w', encoding='utf-8') as handle:
        handle.write('Individual ID\tSex\tAffected\n')
        handle.write('bill\t2\t1\n')
        handle.write('ted\t1\t2\n')

    results = parse_ped_simple(str(tsv))
    assert results == {
        'bill': PedPerson('bill', False, False),
        'ted': PedPerson('ted', True, True),
    }


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
