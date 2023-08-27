"""
tests for the HTML builder
"""
import json
from reanalysis.html_builder import check_date_filter


def test_check_date_filter(tmp_path):
    """

    Returns:

    """
    result_dict = {
        'metadata': {'categories': 'CATEGORY_META'},
        'results': {
            'sample1': {
                'metadata': 'metadata',
                'variants': [{'first_seen': '2021-01-01', 'keep': True}],
            },
            'sample2': {
                'metadata': 'metadata',
                'variants': [{'first_seen': '2022-02-02', 'keep': False}],
            },
        },
    }
    result_path = str(tmp_path / 'results.json')
    with open(result_path, 'w', encoding='utf-8') as handle:
        json.dump(result_dict, handle)

    filtered_results = check_date_filter(result_path, '2021-01-01')
    assert 'sample1' in filtered_results['results']
    assert 'sample2' not in filtered_results['results']
    assert filtered_results['results']['sample1']['variants'][0]['keep']
