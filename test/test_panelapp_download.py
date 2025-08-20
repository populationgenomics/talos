from talos.DownloadPanelApp import (
    PANELS_ENDPOINT,
    get_panels_and_hpo_terms,
    parse_panel,
    parse_panel_activity,
)
from talos.models import HpoTerm


def test_panel_hpo_query(httpx_mock, panels_and_hpos):
    """check that the default parsing delivers correct data"""

    httpx_mock.add_response(url=PANELS_ENDPOINT, json=panels_and_hpos)

    parsed_response = get_panels_and_hpo_terms()

    assert parsed_response == {
        3149: [HpoTerm(id='HP:0011516', label='')],
        4059: [
            HpoTerm(id='HP:0001638', label=''),
            HpoTerm(id='HP:0001637', label=''),
            HpoTerm(id='HP:0011675', label=''),
        ],
        3302: [],
    }


def test_activity_parser(panel_activities):
    """
    check that we correctly parse the activities JSON
    """

    activity_dict = parse_panel_activity(panel_activities)

    # this should be absent
    assert 'NOT_GENE' not in activity_dict

    assert 'GENE1' in activity_dict
    assert activity_dict['GENE1'] == '2022-02-01'

    assert 'GENE2' in activity_dict
    assert activity_dict['GENE2'] == '2024-04-25'

    assert 'GENE3' in activity_dict
    assert activity_dict['GENE3'] == '2023-08-15'


def test_parse_panel(latest_mendeliome, panel_activities):
    result = parse_panel(panel_data=latest_mendeliome, panel_activities=panel_activities)
    assert result == {
        'ENSG00ABCD': {
            'symbol': 'ABCD',
            'chrom': '1',
            'mane_symbol': '',
            'moi': 'biallelic',
            'green_date': '1970-01-01',
        },
        'ENSG00EFGH': {
            'symbol': 'EFGH',
            'chrom': '1',
            'mane_symbol': '',
            'moi': 'monoallelic',
            'green_date': '1970-01-01',
        },
        'ENSG00IJKL': {
            'symbol': 'IJKL',
            'chrom': '1',
            'mane_symbol': '',
            'moi': 'both',
            'green_date': '1970-01-01',
        },
    }


def test_parse_panel_with_mane(latest_mendeliome, panel_activities):
    """check that the default parsing delivers correct data"""
    fake_mane_symbols = {'ABCD': 'EasyAs123D'}
    result = parse_panel(panel_data=latest_mendeliome, panel_activities=panel_activities, symbol_dict=fake_mane_symbols)
    assert result == {
        'EasyAs123D': {
            'symbol': 'ABCD',
            'chrom': '1',
            'mane_symbol': '',
            'moi': 'biallelic',
            'green_date': '1970-01-01',
        },
        'ENSG00ABCD': {
            'symbol': 'ABCD',
            'chrom': '1',
            'mane_symbol': '',
            'moi': 'biallelic',
            'green_date': '1970-01-01',
        },
        'ENSG00EFGH': {
            'symbol': 'EFGH',
            'chrom': '1',
            'mane_symbol': '',
            'moi': 'monoallelic',
            'green_date': '1970-01-01',
        },
        'ENSG00IJKL': {
            'symbol': 'IJKL',
            'chrom': '1',
            'mane_symbol': '',
            'moi': 'both',
            'green_date': '1970-01-01',
        },
    }
