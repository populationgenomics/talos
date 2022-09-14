"""
docstring
"""

import json
import os

from additional_findings.additional_findings_parser import (
    get_data_from_row,
    main,
    triple_protein_changes,
)


def test_triple_protein():
    """
    tests the single to triple IUPAC conversion method
    """

    assert triple_protein_changes(['p.S69R']) == ['p.Ser69Arg']
    assert triple_protein_changes(['p.S69R', 'p.M420L']) == [
        'p.Ser69Arg',
        'p.Met420Leu',
    ]


def test_main_1(tmpdir, acmg_1, acmg_1_exp):
    """
    1st single line test - no specific type
    """
    output = os.path.join(tmpdir, 'output.json')
    main(input_file=acmg_1, output_file=output)
    with open(output, encoding='utf-8') as output_handle:
        actual = json.load(output_handle)
        with open(acmg_1_exp, encoding='utf-8') as expected_handle:
            expected = json.load(expected_handle)
            assert actual == expected, json.dumps(actual, indent=True)


def test_main_2(tmpdir, acmg_2, acmg_2_exp):
    """
    2nd single line test - specific type of change
    """
    output = os.path.join(tmpdir, 'output.json')
    main(input_file=acmg_2, output_file=output)
    with open(output, encoding='utf-8') as output_handle:
        actual = json.load(output_handle)
        with open(acmg_2_exp, encoding='utf-8') as expected_handle:
            expected = json.load(expected_handle)
            assert actual == expected, json.dumps(actual, indent=True)


def test_main_3(tmpdir, acmg_3, acmg_3_exp):
    """
    3rd single line test - specific protein change
    """
    output = os.path.join(tmpdir, 'output.json')
    main(input_file=acmg_3, output_file=output)
    with open(output, encoding='utf-8') as output_handle:
        actual = json.load(output_handle)
        with open(acmg_3_exp, encoding='utf-8') as expected_handle:
            expected = json.load(expected_handle)
            assert actual == expected, json.dumps(actual, indent=True)


def test_main_4(tmpdir, acmg_4, acmg_4_exp):
    """
    4th test, 2 lines, one specific one not
    this test passes, but I don't think it should
    Possible that this type of data won't be provided
    """
    output = os.path.join(tmpdir, 'output.json')
    main(input_file=acmg_4, output_file=output)
    with open(output, encoding='utf-8') as output_handle:
        actual = json.load(output_handle)
        with open(acmg_4_exp, encoding='utf-8') as expected_handle:
            expected = json.load(expected_handle)
            assert actual == expected, json.dumps(actual, indent=True)


def test_data_from_row():
    """
    dictionary test
    """
    data = {
        'Gene': 'ABCDE',
        'Inheritance': 'AR',
        'Disease/Phentyope': 'foo',
        'Phenotype Category': 'bar',
        'Variants to report': 'All P and LP',
    }
    assert get_data_from_row(row_data=data) == {
        'symbol': 'ABCDE',
        'moi': 'Biallelic',
        'flags': ['bar', 'foo'],
        'specific_type': [],
        'specific_variant': [],
    }


def test_data_from_row_2():
    """
    dictionary test
    """
    data = {
        'Gene': 'ABCDE',
        'Inheritance': 'AR',
        'Disease/Phentyope': 'foo',
        'Phenotype Category': 'bar',
        'Variants to report': 'P and LP (pickled variants only)',
    }
    assert get_data_from_row(row_data=data) == {
        'symbol': 'ABCDE',
        'moi': 'Biallelic',
        'flags': ['bar', 'foo'],
        'specific_type': ['pickled'],
        'specific_variant': [],
    }


def test_data_from_row_3():
    """
    dictionary test
    """
    data = {
        'Gene': 'ABCDE',
        'Inheritance': 'AD',
        'Disease/Phentyope': 'foo',
        'Phenotype Category': 'bar',
        'Variants to report': 'p.E69R variants only',
    }
    assert get_data_from_row(row_data=data) == {
        'symbol': 'ABCDE',
        'moi': 'Monoallelic',
        'flags': ['bar', 'foo'],
        'specific_type': [],
        'specific_variant': ['p.Glu69Arg'],
    }
