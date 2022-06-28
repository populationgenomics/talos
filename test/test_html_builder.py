"""
tests for the HTML generation script

Not a huge number, this is not intended as a permanent solution
"""

from reanalysis.html_builder import category_strings, COLORS


def test_numerical_categories():
    """

    :return:
    """
    assert category_strings({}, sample='') == []
    var_data = {'category_2': True}
    assert category_strings(var_data, sample='') == ['2']
    var_data = {'category_4': ['sam1']}
    assert category_strings(var_data, sample='') == []
    assert category_strings(var_data, sample='sam1') == ['de_novo']


def test_set_up_categories():
    """

    :return:
    """
    assert COLORS['de_novo'] == '#FF0000'
