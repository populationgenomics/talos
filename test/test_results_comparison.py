"""
test file for annotation with first-seen dates
"""

from copy import deepcopy
from datetime import datetime
from os.path import join
from time import sleep

import pytest
import zoneinfo

from reanalysis.models import (
    Coordinates,
    HistoricSampleVariant,
    HistoricVariants,
    ReportVariant,
    ResultData,
    SmallVariant,
)
from reanalysis.static_values import get_granular_date
from reanalysis.utils import date_annotate_results, date_from_string, find_latest_file

COORD_1 = Coordinates(chrom='1', pos=1, ref='A', alt='G')
COORD_2 = Coordinates(chrom='2', pos=2, ref='A', alt='G')

VAR_1 = SmallVariant(coordinates=COORD_1, info={}, transcript_consequences=[])
VAR_2 = SmallVariant(coordinates=COORD_2, info={}, transcript_consequences=[])
REP_SAM1_1 = ReportVariant(sample='sam1', var_data=VAR_1, categories={'1'}, gene='ENSG1')
REP_SAM1_1_Independent = ReportVariant(sample='sam1', var_data=VAR_1, categories={'1'}, gene='ENSG1', independent=True)
REP_SAM2_2 = ReportVariant(sample='sam2', var_data=VAR_2, categories={'2'}, gene='ENSG2', independent=True)
REP_SAM2_12 = ReportVariant(sample='sam2', var_data=VAR_2, categories={'1', '2'}, gene='ENSG2')

TIMEZONE = zoneinfo.ZoneInfo('Australia/Brisbane')
OLD_DATE = datetime(year=2000, month=1, day=1, tzinfo=TIMEZONE).strftime('%Y-%m-%d')


def test_date_annotate_one():
    """
    same category, old date - revert 'first-seen'
    """
    results = ResultData(
        **{
            'results': {
                'sam1': {
                    'variants': [deepcopy(REP_SAM1_1)],
                    'metadata': {'ext_id': 'jeff', 'family_id': 'jeff'},
                },
            },
        },
    )
    historic = HistoricVariants(
        **{
            'results': {
                'sam1': {
                    VAR_1.coordinates.string_format: {
                        'categories': {'1': OLD_DATE},
                        'independent': False,
                        'first_tagged': OLD_DATE,
                    },
                },
            },
        },
    )
    date_annotate_results(results, historic)
    assert results.results['sam1'].variants[0].first_tagged == OLD_DATE


def test_date_annotate_two():
    """
    if there's a new category, don't update dates
    """
    results = ResultData(
        **{
            'results': {
                'sam2': {
                    'metadata': {'ext_id': 'sam2', 'family_id': '2'},
                    'variants': [deepcopy(REP_SAM2_12)],
                },
            },
        },
    )
    prior_results = deepcopy(results)
    historic = HistoricVariants(
        **{
            'results': {
                'sam2': {
                    COORD_2.string_format: {
                        'categories': {'1': OLD_DATE},
                        'independent': False,
                        'first_tagged': get_granular_date(),
                    },
                },
            },
        },
    )
    date_annotate_results(results, historic)
    assert len(historic.results['sam2'][COORD_2.string_format].categories) == 2
    assert historic.results['sam2'][COORD_2.string_format].categories['1'] == OLD_DATE
    assert historic.results['sam2'][COORD_2.string_format].categories['2'] == get_granular_date()
    assert results == prior_results


def test_date_annotate_three():
    """
    if a variant is newly independent, DON'T update the dates
    """
    results = ResultData(
        **{
            'results': {
                'sam1': {
                    'metadata': {'ext_id': 'sam1', 'family_id': '1'},
                    'variants': [deepcopy(REP_SAM1_1_Independent)],
                },
            },
        },
    )
    historic = HistoricVariants(
        **{
            'results': {
                'sam1': {
                    COORD_1.string_format: {
                        'categories': {'1': OLD_DATE},
                        'support_vars': ['flipflop'],
                        'independent': False,
                        'first_tagged': OLD_DATE,
                    },
                },
            },
        },
    )
    date_annotate_results(results, historic)
    assert historic.results['sam1'][COORD_1.string_format] == HistoricSampleVariant(
        **{
            'categories': {'1': OLD_DATE},
            'support_vars': ['flipflop'],
            'independent': True,
            'first_tagged': OLD_DATE,
        },
    )


def test_date_annotate_four():
    """
    if there's a new category, don't update dates
    """
    results = ResultData(
        **{
            'results': {
                'sam1': {
                    'metadata': {'ext_id': 'sam1', 'family_id': '1'},
                    'variants': [deepcopy(REP_SAM1_1)],
                },
                'sam2': {
                    'metadata': {'ext_id': 'sam2', 'family_id': '2'},
                    'variants': [deepcopy(REP_SAM2_2)],
                },
            },
        },
    )
    historic = HistoricVariants(
        **{
            'results': {
                'sam1': {COORD_1.string_format: {'categories': {'2': OLD_DATE}, 'first_tagged': get_granular_date()}},
                'sam3': {},
            },
        },
    )
    date_annotate_results(results, historic)
    assert historic == HistoricVariants(
        **{
            'results': {
                'sam1': {
                    COORD_1.string_format: {
                        'categories': {'1': get_granular_date(), '2': OLD_DATE},
                        'first_tagged': get_granular_date(),
                    },
                },
                'sam2': {
                    COORD_2.string_format: {
                        'categories': {'2': get_granular_date()},
                        'first_tagged': get_granular_date(),
                    },
                },
                'sam3': {},
            },
        },
    )


def touch(filepath: str):
    """simulates a touch"""
    with open(filepath, 'w', encoding='utf-8'):
        pass


def test_find_latest(tmp_path):
    """
    check that we find the correct latest file
    """
    tmp_str = str(tmp_path)

    touch(join(tmp_str, '2025-10-10.json'))
    sleep(0.2)
    touch(join(tmp_str, '2026-10-10.json'))
    sleep(0.2)
    touch(join(tmp_str, '2024-10-10.json'))
    assert '2026-10-10.json' in find_latest_file(results_folder=tmp_str, dataset='cohort')


def test_find_latest_singletons(tmp_path):
    """
    check that we find the correct latest file
    """
    tmp_str = str(tmp_path)
    touch(join(tmp_str, 'singletons_2020-10-10.json'))
    sleep(0.2)
    touch(join(tmp_str, '2025-10-10.json'))
    sleep(0.2)
    touch(join(tmp_str, '2025-11-10.json'))
    assert 'singletons_2020-10-10.json' in find_latest_file(
        results_folder=tmp_str,
        start='singletons',
        dataset='cohort',
    )


def test_find_latest_with_ext(tmp_path):
    """
    check that we find the correct latest file
    """
    tmp_str = str(tmp_path)
    touch(join(tmp_str, '2020-10-10.txt'))
    sleep(0.2)
    touch(join(tmp_str, '2026-10-10.txt'))
    sleep(0.2)
    touch(join(tmp_str, '2027-10-10.json'))
    assert '2026-10-10.txt' in find_latest_file(results_folder=tmp_str, ext='txt', dataset='cohort')


def test_date_recovery():
    """
    check that we can correctly pull the dates from a String
    """

    assert date_from_string('2020-10-10') == '2020-10-10'
    assert date_from_string('2020-10-10T10:10:10') == '2020-10-10'
    assert date_from_string('PICKLED2020-10-10T10:10:10.000') == '2020-10-10'
    with pytest.raises(ValueError):
        date_from_string('PICKLEDNO DATET10:10:10.000')
