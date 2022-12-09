"""
test file for the results subtraction methods in utils
"""

from dataclasses import dataclass, field
from os.path import join
from time import sleep

from reanalysis.utils import (
    subtract_results,
    add_results,
    find_latest,
    Coordinates,
    TODAY,
)


# pylint: disable=consider-using-with


@dataclass
class MiniVariant:
    """
    tiny class to simulate the variants
    """

    categories: list
    coords: Coordinates


@dataclass
class MiniReport:
    """
    tiny class to simulate the reports
    """

    var_data: MiniVariant
    support_vars: list[str] = field(default_factory=list)


COORD_1 = Coordinates('1', 1, 'A', 'G')
COORD_2 = Coordinates('2', 2, 'A', 'G')
COORD_3 = Coordinates('3', 3, 'A', 'G')
GENERIC_REPORT = MiniReport(MiniVariant(categories=['1'], coords=COORD_1))
GENERIC_REPORT_2 = MiniReport(MiniVariant(categories=['2'], coords=COORD_2))
GENERIC_REPORT_3 = MiniReport(MiniVariant(categories=['3'], coords=COORD_3))
GENERIC_REPORT_12 = MiniReport(MiniVariant(categories=['1', '2'], coords=COORD_1))


def test_subtraction_null():
    """
    subtraction if nothing to subtract
    """
    new = {'sample': [GENERIC_REPORT]}
    old = {}
    assert subtract_results(new, old) == new


def test_subtraction_no_matches():
    """
    no previous results for this sample
    """
    new = {'sample': [GENERIC_REPORT]}
    old = {'sample2': [1]}
    assert subtract_results(new, old) == new


def test_subtraction_one_exact():
    """
    one variant fully removed
    """
    new = {
        'sample': [GENERIC_REPORT],
        'sample2': [GENERIC_REPORT],
    }
    old = {
        'sample': {COORD_1.string_format: {'categories': {'1': 1}, 'support_vars': []}}
    }
    assert subtract_results(new, old) == {'sample': [], 'sample2': [GENERIC_REPORT]}


def test_subtraction_match_no_categories():
    """
    variant matches, but categories do not
    """
    new = {'sample': [GENERIC_REPORT]}
    old = {
        'sample': {COORD_1.string_format: {'categories': {'2': 1}, 'support_vars': []}}
    }
    assert subtract_results(new, old) == new


def test_subtraction_partial_categories():
    """
    variant matches, but categories only partially match
    """
    new = {'sample': [MiniReport(MiniVariant(categories=['1', '2'], coords=COORD_1))]}
    old = {'sample': {COORD_1.string_format: {'categories': {'2': 1}}}}
    assert subtract_results(new, old) == {
        'sample': [MiniReport(MiniVariant(categories=['1'], coords=COORD_1))]
    }


def test_subtraction_new_partner():
    """
    variant matches, categories match, but new comp-het
    """
    new = {
        'sample': [
            MiniReport(
                MiniVariant(categories=['1'], coords=COORD_1), support_vars=['foo']
            )
        ]
    }
    old = {
        'sample': {
            COORD_1.string_format: {'categories': {'1': 1}, 'support_vars': ['bar']}
        }
    }
    assert subtract_results(new, old) == {
        'sample': [
            MiniReport(
                MiniVariant(categories=['1'], coords=COORD_1), support_vars=['foo']
            )
        ]
    }


def test_add_new_sample():
    """
    add a novel sample
    """
    cum = {}
    new = {'sample': [GENERIC_REPORT]}
    add_results(new, cum)
    assert cum == {
        'sample': {
            COORD_1.string_format: {'categories': {'1': TODAY}, 'support_vars': []}
        }
    }


def test_add_new_sample_variant():
    """
    integrate a new sample for an existing sample
    """
    new = {'sample': [GENERIC_REPORT]}
    cum = {'sample': {'1': {'2'}}}
    add_results(new, cum)
    assert cum == {
        'sample': {
            COORD_1.string_format: {'categories': {'1': TODAY}, 'support_vars': []},
            '1': {'2'},
        }
    }


def test_add_new_category_vardup():
    """
    integrate a new sample for an existing sample
    """
    # ludicrous situation, but should be manageable
    new = {'sample': [GENERIC_REPORT, GENERIC_REPORT_12]}
    cum = {'sample': {}}
    add_results(new, cum)
    assert cum == {
        'sample': {
            COORD_1.string_format: {
                'categories': {'1': TODAY, '2': TODAY},
                'support_vars': [],
            }
        }
    }


def test_add_support_vars():
    """
    integrate a new sample for an existing sample
    """
    # ludicrous situation, but should be manageable
    new = {
        'sample': [
            GENERIC_REPORT,
            MiniReport(
                MiniVariant(categories=['999'], coords=COORD_1), support_vars=['foobar']
            ),
        ]
    }
    cum = {
        'sample': {
            COORD_1.string_format: {
                'categories': {'1': 1, '2': 1},
                'support_vars': [],
            }
        }
    }
    add_results(new, cum)
    assert cum == {
        'sample': {
            COORD_1.string_format: {
                'categories': {'1': 1, '2': 1, '999': TODAY},
                'support_vars': ['foobar'],
            }
        }
    }


def test_add_various():
    """
    add a bunch of things
    """
    new = {
        'sample': [
            MiniReport(
                MiniVariant(categories=['999'], coords=COORD_1), support_vars=['foobar']
            ),
            MiniReport(
                MiniVariant(categories=['A', 'B'], coords=COORD_1),
                support_vars=[],
            ),
        ],
        'sample2': [GENERIC_REPORT],
        'sample3': [
            MiniReport(
                MiniVariant(categories=['B'], coords=COORD_3), support_vars=['foobar']
            )
        ],
    }
    cum = {
        'sample': {
            COORD_1.string_format: {'categories': {'A': 1, 'C': 2}, 'support_vars': []},
            COORD_2.string_format: {'categories': {'1': 1}, 'support_vars': ['foo']},
            COORD_3.string_format: {
                'categories': {'A': 1, 'B': 2, 'C': 3},
                'support_vars': [],
            },
        },
        'sample3': {
            COORD_3.string_format: {'categories': {'B': 1}, 'support_vars': ['batman']}
        },
    }
    add_results(new, cum)
    assert cum == {
        'sample': {
            COORD_1.string_format: {
                'categories': {'A': 1, 'B': TODAY, 'C': 2, '999': TODAY},
                'support_vars': ['foobar'],
            },
            COORD_2.string_format: {'categories': {'1': 1}, 'support_vars': ['foo']},
            COORD_3.string_format: {
                'categories': {'A': 1, 'B': 2, 'C': 3},
                'support_vars': [],
            },
        },
        'sample2': {
            COORD_1.string_format: {'categories': {'1': TODAY}, 'support_vars': []}
        },
        'sample3': {
            COORD_3.string_format: {
                'categories': {'B': 1},
                'support_vars': ['batman', 'foobar'],
            }
        },
    }


def test_find_latest(tmp_path):
    """
    check that we find the correct latest file
    """
    open(join(str(tmp_path), 'file1.json'), 'w', encoding='utf-8')
    sleep(0.2)
    open(join(str(tmp_path), 'file2.json'), 'w', encoding='utf-8')
    sleep(0.2)
    open(join(str(tmp_path), 'file3.json'), 'w', encoding='utf-8')
    assert 'file3.json' in find_latest(str(tmp_path), singletons=False)


def test_find_latest_singletons(tmp_path):
    """
    check that we find the correct latest file
    """
    open(join(str(tmp_path), 'singletons_file1.json'), 'w', encoding='utf-8')
    sleep(0.2)
    open(join(str(tmp_path), 'file2.json'), 'w', encoding='utf-8')
    sleep(0.2)
    open(join(str(tmp_path), 'file3.json'), 'w', encoding='utf-8')
    assert 'singletons_file1.json' in find_latest(str(tmp_path), singletons=True)


def test_find_latest_with_ext(tmp_path):
    """
    check that we find the correct latest file
    """
    open(join(str(tmp_path), 'file1.txt'), 'w', encoding='utf-8')
    sleep(0.2)
    open(join(str(tmp_path), 'file2.txt'), 'w', encoding='utf-8')
    sleep(0.2)
    open(join(str(tmp_path), 'file3.txt'), 'w', encoding='utf-8')
    assert 'file3.txt' in find_latest(str(tmp_path), ext='txt', singletons=False)
