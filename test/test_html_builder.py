"""
tests for the HTML builder
"""

from test.test_utils import TWO_EXPECTED

from talos.original_templates.CreateTalosHTML import check_date_filter
from talos.models import Coordinates, ResultData, SmallVariant

TEST_COORDS = Coordinates(chrom='1', pos=1, ref='A', alt='C')
VAR_1 = SmallVariant(coordinates=TEST_COORDS, info={}, transcript_consequences=[])
TEST_COORDS_2 = Coordinates(chrom='2', pos=2, ref='A', alt='C')
VAR_2 = SmallVariant(coordinates=TEST_COORDS_2, info={}, transcript_consequences=[])


def test_check_date_filter(tmp_path):
    """

    Returns:

    """
    # needs to be a valid ResultData
    result_dict = ResultData(
        metadata={'categories': {'1': '1'}, 'run_datetime': '2021-01-01'},
        results={
            'sample1': {
                'metadata': {'ext_id': 'sample1', 'family_id': 'sample1'},
                'variants': [
                    {
                        'sample': 'sample1',
                        'first_tagged': '2021-01-01',
                        'var_data': VAR_1,
                    },
                ],
            },
            'sample2': {
                'metadata': {'ext_id': 'sample1', 'family_id': 'sample1'},
                'variants': [
                    {
                        'sample': 'sample1',
                        'first_tagged': '2022-02-02',
                        'var_data': VAR_2,
                    },
                ],
            },
        },
    )
    result_path = str(tmp_path / 'results.json')
    with open(result_path, 'w', encoding='utf-8') as handle:
        handle.write(ResultData.model_validate(result_dict).model_dump_json(indent=4))

    filtered_results = check_date_filter(result_path, '2021-01-01')
    assert 'sample1' in filtered_results.results
    assert 'sample2' not in filtered_results.results
    assert filtered_results.results['sample1'].variants[0].first_tagged == '2021-01-01'
    assert filtered_results.results['sample1'].variants[0].var_data.coordinates.string_format == '1-1-A-C'


def test_check_date_filter_pair(tmp_path):
    """
    test for the case when a comp-het pair is both recovered by one meeting the
    specified date and the other not
    """
    # needs to be a valid ResultData
    result_dict = ResultData(
        metadata={'categories': {'1': '1'}, 'run_datetime': '2021-01-01'},
        results={
            'sample1': {
                'metadata': {'ext_id': 'sample1', 'family_id': 'sample1'},
                'variants': [
                    {
                        'sample': 'sample1',
                        'first_tagged': '2021-01-01',
                        'var_data': VAR_1,
                        'support_vars': {VAR_2.coordinates.string_format},
                    },
                    {
                        'sample': 'sample1',
                        'first_tagged': '2023-03-03',
                        'var_data': VAR_2,
                        'support_vars': {VAR_1.coordinates.string_format},
                    },
                ],
            },
            'sample2': {
                'metadata': {'ext_id': 'sample1', 'family_id': 'sample1'},
                'variants': [
                    {
                        'sample': 'sample1',
                        'first_tagged': '2022-02-02',
                        'var_data': VAR_2,
                    },
                ],
            },
        },
    )
    result_path = str(tmp_path / 'results.json')
    with open(result_path, 'w', encoding='utf-8') as handle:
        handle.write(ResultData.model_validate(result_dict).model_dump_json(indent=4))

    filtered_results = check_date_filter(result_path, '2021-01-01')
    assert 'sample1' in filtered_results.results
    assert 'sample2' not in filtered_results.results
    assert len(filtered_results.results['sample1'].variants) == TWO_EXPECTED
    assert sorted(var.first_tagged for var in filtered_results.results['sample1'].variants) == [
        '2021-01-01',
        '2023-03-03',
    ]


def test_check_date_filter_none_pass(tmp_path):
    """ """
    # needs to be a valid ResultData
    result_dict = ResultData(
        metadata={'categories': {'1': '1'}, 'run_datetime': '2023-01-01'},
        results={
            'sample1': {
                'metadata': {'ext_id': 'sample1', 'family_id': 'sample1'},
                'variants': [
                    {
                        'sample': 'sample1',
                        'first_tagged': '2021-01-01',
                        'var_data': VAR_1,
                    },
                ],
            },
            'sample2': {
                'metadata': {'ext_id': 'sample2', 'family_id': 'sample2'},
                'variants': [
                    {
                        'sample': 'sample2',
                        'first_tagged': '2022-02-02',
                        'var_data': VAR_2,
                    },
                ],
            },
        },
    )
    result_path = str(tmp_path / 'results.json')
    with open(result_path, 'w', encoding='utf-8') as handle:
        handle.write(ResultData.model_validate(result_dict).model_dump_json(indent=4))

    assert check_date_filter(result_path, '2023-01-01') is None


def test_check_date_filter_date_from_meta(tmp_path):
    """ """
    # needs to be a valid ResultData
    result_dict = ResultData(
        metadata={'categories': {'1': '1'}, 'run_datetime': '2022-02-02'},
        results={
            'sample1': {
                'metadata': {'ext_id': 'sample1', 'family_id': 'sample1'},
                'variants': [
                    {
                        'sample': 'sample1',
                        'first_tagged': '2021-01-01',
                        'var_data': VAR_1,
                    },
                ],
            },
            'sample2': {
                'metadata': {'ext_id': 'sample1', 'family_id': 'sample1'},
                'variants': [
                    {
                        'sample': 'sample2',
                        'first_tagged': '2022-02-02',
                        'var_data': VAR_2,
                    },
                ],
            },
        },
    )
    result_path = str(tmp_path / 'results.json')
    with open(result_path, 'w', encoding='utf-8') as handle:
        handle.write(ResultData.model_validate(result_dict).model_dump_json(indent=4))

    filtered_results = check_date_filter(result_path)
    assert 'sample1' not in filtered_results.results
    assert 'sample2' in filtered_results.results
    assert filtered_results.results['sample2'].variants[0].first_tagged == '2022-02-02'
    assert filtered_results.results['sample2'].variants[0].var_data.coordinates.string_format == '2-2-A-C'
