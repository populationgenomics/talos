"""
tests for the HTML generation script

Not a huge number, this is not intended as a permanent solution
"""

from reanalysis.html_builder import COLORS


def test_set_up_categories():
    """

    :return:
    """
    assert COLORS['de_novo'] == '#FF0000'
