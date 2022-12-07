"""
tests for the PanelApp parser
"""

import json
import os
from unittest.mock import patch

import pytest

from reanalysis.new_panelapp_query import (
    get_panel_green,
    is_this_gene_new,
    read_panels_from_participant_file,
)


PWD = os.path.dirname(__file__)
INPUT = os.path.join(PWD, 'input')
PANELAPP_LATEST = os.path.join(INPUT, 'panelapp_current_137.json')
PANELAPP_INCIDENTALOME = os.path.join(INPUT, 'incidentalome.json')
PANELAPP_OLDER = os.path.join(INPUT, 'panelapp_older_137.json')
LATEST_EXPECTED = os.path.join(INPUT, 'panel_green_latest_expected.json')
CHANGES_EXPECTED = os.path.join(INPUT, 'panel_changes_expected.json')
FAKE_GENE_LIST = os.path.join(INPUT, 'fake_gene_list.json')
OLD_VERSION = 'old'


@pytest.fixture(name='fake_panelapp')
def fixture_fake_panelapp(requests_mock):
    """
    a new fixture to contain the panel data
    :param requests_mock:
    :return:
    """
    with open(PANELAPP_LATEST, 'r', encoding='utf-8') as handle:
        requests_mock.register_uri(
            'GET',
            'https://panelapp.agha.umccr.org/api/v1/panels/137',
            json=json.load(handle),
        )
    with open(PANELAPP_INCIDENTALOME, 'r', encoding='utf-8') as handle:
        requests_mock.register_uri(
            'GET',
            'https://panelapp.agha.umccr.org/api/v1/panels/126',
            json=json.load(handle),
        )


@patch(
    'reanalysis.new_panelapp_query.OLD_DATA',
    {'genes': {'ENSG00ABCD': {'panels': ['all']}, 'ENSG00EFGH': {'panels': ['all']}}},
)
def test_panel_query(fake_panelapp):  # pylint: disable=unused-argument
    """
    check that the default parsing delivers correct data
    oof, this was a tricky one
    :param fake_panelapp: fake web hook mock
    """

    gd = {'genes': {}, 'metadata': []}
    get_panel_green(gd, panel_id=137)
    with open(LATEST_EXPECTED, 'r', encoding='utf-8') as handle:
        expected = json.load(handle)

    print(gd)
    print(expected)
    assert gd == expected


@patch(
    'reanalysis.new_panelapp_query.OLD_DATA',
    {
        'genes': {
            'ENSG00ABCD': {'panels': [137]},
            'ENSG00EFGH': {'panels': [137]},
            'ENSG00IJKL': {'panels': [137]},
        },
    },
)
def test_panel_query_addition(
    fake_panelapp, panel_updates
):  # pylint: disable=unused-argument
    """
    check that the default parsing delivers correct data
    oof, this was a tricky one
    :param fake_panelapp: fake web hook mock
    """

    gd = {
        'metadata': [{'version': '0.11088', 'name': 'Mendeliome', 'id': 137}],
        'genes': {
            'ENSG00ABCD': {
                'symbol': 'ABCD',
                'moi': 'Monoallelic',
                'new': True,
                'panels': [137],
            },
            'ENSG00IJKL': {
                'symbol': 'IJKL',
                'moi': 'Mono_And_Biallelic',
                'new': True,
                'panels': [137],
            },
        },
    }

    get_panel_green(gd, panel_id=126)
    assert gd == panel_updates


#
#
# def test_gene_list_changes():
#     """
#     check usage of a gene list
#     :return:
#     """
#
#     with open(LATEST_EXPECTED, 'r', encoding='utf-8') as handle:
#         latest = json.load(handle)
#     previous_genes = {'ABCD'}
#     gene_list_differences(latest, previous_genes)
#     assert not latest['ENSG00ABCD']['new']
#     assert latest['ENSG00IJKL']['new']
#
#
# def test_get_json_response(fake_panelapp):  # pylint: disable=unused-argument
#     """
#     read the json content via an API call
#     :param fake_panelapp:
#     :return:
#     """
#     result = get_json_response('https://panelapp.agha.umccr.org/api/v1/panels/137')
#     with open(PANELAPP_LATEST, 'r', encoding='utf-8') as handle:
#         assert result == json.load(handle)
#
#
# def test_parse_local_gene_list():
#     """
#     test parsing of a local file
#     :return:
#     """
#     found_genes = parse_gene_list(FAKE_GENE_LIST)
#     assert found_genes == {'foo bar', 'foo', 'bar'}
#
#
# def test_second_panel_update(fake_panelapp):  # pylint: disable=unused-argument
#     """
#     check that the combination method works
#     """
#     main = get_panel_green('137')
#     early_keys = set(main.keys())
#     additional = get_panel_green('126')
#     combine_mendeliome_with_other_panels(main, additional)
#     assert set(main.keys()) == early_keys
#     assert main['ENSG00ABCD'].get('flags') == ['Incidentalome']
#
#
# def test_second_panel_update_moi():
#     """
#     panels overlap with MOI in second panel
#     - expect overwriting of MOI
#     """
#     main = {
#         'metadata': [{'name': 'NAME', 'id': 'fish'}],
#         'ensg1': {'entity_name': '1', 'moi': None, 'flags': []},
#     }
#     additional = {
#         'metadata': [{'name': 'NAME', 'version': '1.0', 'id': '123'}],
#         'ensg1': {'entity_name': '1', 'moi': 'REALLY_BIG'},
#     }
#     combine_mendeliome_with_other_panels(main, additional)
#     assert main['ensg1'].get('flags') == ['NAME']
#     assert main['ensg1'].get('moi') == 'REALLY_BIG'
#     assert len(main['metadata']) == 2
#     assert main['metadata'][0]['id'] == 'fish'
#     assert main['metadata'][1]['id'] == '123'


def test_get_list_from_participants(tmp_path):
    """
    tests the unique panel finder
    """
    party_data = {
        'i': {'panels': [1, 2], 'what': 'does'},
        'am': {'panels': [1, 3], 'the': 'fox'},
        'sam': {'panels': [9, 99], 'say?': 'Wa-pa-pa-pa-pa-pa-pow!'},
    }
    tmp_json = tmp_path / 'temp.json'
    with open(tmp_json, 'w', encoding='utf-8') as handle:
        json.dump(party_data, handle)
    assert read_panels_from_participant_file(tmp_json) == {1, 2, 3, 9, 99}


@patch('reanalysis.new_panelapp_query.OLD_DATA', {'genes': {'ABC': {'panels': [1]}}})
def test_do_stuff():
    """
    mock in the historical data and check for when a gene is new
    """
    assert is_this_gene_new('ABC', 12)
    assert not is_this_gene_new('ABC', 1)


@patch(
    'reanalysis.new_panelapp_query.OLD_DATA', {'genes': {'ABC': {'panels': ['all']}}}
)
def test_do_stuff_all():
    """
    historical data says this gene is perceived as in 'all panels'
    """
    assert not is_this_gene_new('ABC', 12)
    assert not is_this_gene_new('ABC', 1)
