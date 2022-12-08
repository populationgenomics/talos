"""
tests for the PanelApp parser
"""

import json
from unittest.mock import patch

import pytest

from reanalysis.query_panelapp import (
    get_panel_green,
    is_this_gene_new,
    read_panels_from_participant_file,
)


@pytest.fixture(name='fake_panelapp')
def fixture_fake_panelapp(requests_mock, latest_mendeliome, latest_incidentalome):
    """
    prepares the web requests mock to serve as stand-in panelapp
    Args:
        requests_mock ():
        latest_mendeliome ():
        latest_incidentalome ():
    """

    requests_mock.register_uri(
        'GET',
        'https://panelapp.agha.umccr.org/api/v1/panels/137',
        json=latest_mendeliome,
    )
    requests_mock.register_uri(
        'GET',
        'https://panelapp.agha.umccr.org/api/v1/panels/126',
        json=latest_incidentalome,
    )


@patch(
    'reanalysis.new_panelapp_query.OLD_DATA',
    {'genes': {'ENSG00ABCD': {'panels': ['all']}, 'ENSG00EFGH': {'panels': ['all']}}},
)
def test_panel_query(
    fake_panelapp, mendeliome_expected
):  # pylint: disable=unused-argument
    """
    check that the default parsing delivers correct data
    oof, this was a tricky one
    :param fake_panelapp: fake web hook mock
    """

    gd = {'genes': {}, 'metadata': []}
    get_panel_green(gd, panel_id=137)
    assert gd == mendeliome_expected


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
    # assumed data we already gathered
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

    # should query for and integrate the incidentalome content
    get_panel_green(gd, panel_id=126)
    assert gd == panel_updates


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
