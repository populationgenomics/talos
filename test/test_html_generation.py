"""
tests for the results -> HTML module
"""
from typing import List
import pytest

from reanalysis.html_builder import numerical_categories, string_format_coords


@pytest.mark.parametrize(
    'one,two,three,four,result',
    (
        (True, True, True, True, ['1', '2', '3', '4']),
        (True, False, False, False, ['1']),
        (False, True, False, False, ['2']),
        (False, False, True, True, ['3', '4']),
        (False, False, False, False, []),
    ),
)
def test_numerical_categories(
    one: bool, two: bool, three: bool, four: bool, result: List[str]
):
    """
    ronseal
    :return:
    """

    var_data = {
        'category_1': one,
        'category_2': two,
        'category_3': three,
        'category_4': four,
    }
    assert numerical_categories(var_data) == result


def test_string_format_coords():
    """
    ronseal
    :return:
    """
    coords = {
        'chrom': 'N',
        'pos': 'I',
        'ref': 'C',
        'alt': 'E',
    }
    assert string_format_coords(coords) == 'N-I-C-E'
