"""
tests for the HTML generation script

Not a huge number, this is not intended as a permanent solution
"""

from reanalysis.html_builder import numerical_categories, set_up_colors


def test_numerical_categories():
    """

    :return:
    """
    assert numerical_categories({}, sample='') == []
    var_data = {'category_2': True}
    assert numerical_categories(var_data, sample='') == ['2']
    var_data = {'category_4': ['sam1']}
    assert numerical_categories(var_data, sample='') == []
    assert numerical_categories(var_data, sample='sam1') == ['de_novo']


def test_set_up_categories():
    """

    :return:
    """
    cols = set_up_colors()
    assert cols['de_novo'] == '#FF0000'
