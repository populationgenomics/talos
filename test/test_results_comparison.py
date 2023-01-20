"""
test file for annotation with first-seen dates
"""


from dataclasses import dataclass, field
from datetime import datetime
from os.path import join
from time import sleep
from copy import deepcopy

from reanalysis.utils import (
    date_annotate_results,
    find_latest_file,
    Coordinates,
    GRANULAR,
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
    first_seen: str = GRANULAR


COORD_1 = Coordinates('1', 1, 'A', 'G')
COORD_2 = Coordinates('2', 2, 'A', 'G')

GENERIC_REPORT = MiniReport(MiniVariant(categories=['1'], coords=COORD_1))
GENERIC_REPORT_12 = MiniReport(MiniVariant(categories=['1', '2'], coords=COORD_1))
GENERIC_REPORT_2 = MiniReport(MiniVariant(categories=['2'], coords=COORD_2))

OLD_DATE = datetime(year=2000, month=1, day=1).strftime('%Y-%m-%d')


def test_date_annotate_one():
    """
    checks we can add one entry to a None historic
    """
    results = {'sample': {'variants': [GENERIC_REPORT]}}
    new_results, cumulative = date_annotate_results(results)
    assert results == new_results
    assert cumulative == {
        'sample': {
            COORD_1.string_format: {'categories': {'1': GRANULAR}, 'support_vars': []}
        }
    }


def test_date_annotate_two():
    """
    same category, old date - revert 'first-seen'
    """
    results = {'sample': {'variants': [deepcopy(GENERIC_REPORT)]}}
    historic = {
        'sample': {
            COORD_1.string_format: {'categories': {'1': OLD_DATE}, 'support_vars': []}
        }
    }
    results, cumulative = date_annotate_results(results, historic)
    assert cumulative == historic
    assert results == {
        'sample': {
            'variants': [
                MiniReport(
                    MiniVariant(categories=['1'], coords=COORD_1), first_seen=OLD_DATE
                )
            ]
        }
    }


def test_date_annotate_three():
    """
    if there's a new category, don't update dates
    """
    results = {'sample': {'variants': [deepcopy(GENERIC_REPORT_12)]}}
    historic = {
        'sample': {
            COORD_1.string_format: {'categories': {'1': OLD_DATE}, 'support_vars': []}
        }
    }
    new_results, cumulative = date_annotate_results(results, historic)
    assert cumulative == {
        'sample': {
            COORD_1.string_format: {
                'categories': {'1': OLD_DATE, '2': GRANULAR},
                'support_vars': [],
            }
        }
    }
    assert results == new_results


def test_date_annotate_four():
    """
    if there's a new category, don't update dates
    """
    results = {'sample': {'variants': [deepcopy(GENERIC_REPORT)]}}
    historic = {
        'sample': {
            COORD_1.string_format: {
                'categories': {'1': OLD_DATE},
                'support_vars': ['flipflop'],
            }
        }
    }
    new_results, historic = date_annotate_results(results, historic)
    assert historic == {
        'sample': {
            COORD_1.string_format: {
                'categories': {'1': OLD_DATE},
                'support_vars': ['flipflop'],
            }
        }
    }
    assert results == new_results


def test_date_annotate_five():
    """
    if there's a new category, don't update dates
    """
    results = {
        'sample': {'variants': [deepcopy(GENERIC_REPORT)]},
        'sample2': {'variants': [deepcopy(GENERIC_REPORT_2)]},
    }
    historic = {
        'sample': {
            COORD_1.string_format: {
                'categories': {'2': OLD_DATE},
                'support_vars': [],
            }
        },
        'sample3': {},
    }
    results, historic = date_annotate_results(results, historic)
    assert historic == {
        'sample': {
            COORD_1.string_format: {
                'categories': {'1': GRANULAR, '2': OLD_DATE},
                'support_vars': [],
            }
        },
        'sample2': {
            COORD_2.string_format: {
                'categories': {'2': GRANULAR},
                'support_vars': [],
            }
        },
        'sample3': {},
    }
    assert results == {
        'sample': {
            'variants': [
                MiniReport(
                    MiniVariant(categories=['1'], coords=COORD_1),
                    first_seen=GRANULAR,
                )
            ]
        },
        'sample2': {
            'variants': [
                MiniReport(
                    MiniVariant(categories=['2'], coords=COORD_2),
                    first_seen=GRANULAR,
                )
            ]
        },
    }


def test_find_latest(tmp_path):
    """
    check that we find the correct latest file
    """
    tmp_str = str(tmp_path)

    open(join(tmp_str, 'file1.json'), 'w', encoding='utf-8')
    sleep(0.2)
    open(join(tmp_str, 'file2.json'), 'w', encoding='utf-8')
    sleep(0.2)
    open(join(tmp_str, 'file3.json'), 'w', encoding='utf-8')
    assert 'file3.json' in find_latest_file(tmp_str)


def test_find_latest_singletons(tmp_path):
    """
    check that we find the correct latest file
    """
    tmp_str = str(tmp_path)
    open(join(tmp_str, 'singletons_file1.json'), 'w', encoding='utf-8')
    sleep(0.2)
    open(join(tmp_str, 'file2.json'), 'w', encoding='utf-8')
    sleep(0.2)
    open(join(tmp_str, 'file3.json'), 'w', encoding='utf-8')
    assert 'singletons_file1.json' in find_latest_file(tmp_str, start='singletons')


def test_find_latest_with_ext(tmp_path):
    """
    check that we find the correct latest file
    """
    tmp_str = str(tmp_path)
    open(join(tmp_str, 'file1.txt'), 'w', encoding='utf-8')
    sleep(0.2)
    open(join(tmp_str, 'file2.txt'), 'w', encoding='utf-8')
    sleep(0.2)
    open(join(tmp_str, 'file3.txt'), 'w', encoding='utf-8')
    assert 'file3.txt' in find_latest_file(tmp_str, ext='txt')
