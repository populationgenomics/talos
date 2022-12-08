"""
tests for the PanelApp parser
"""

import json
from unittest.mock import patch

import pytest

from reanalysis.query_panelapp import (
    get_panel_green,
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
    'reanalysis.query_panelapp.OLD_DATA',
    {'genes': {'ENSG00ABCD': {'panels': ['pink']}, 'ENSG00EFGH': {'panels': [137]}}},
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
    'reanalysis.query_panelapp.OLD_DATA',
    {
        'genes': {
            'ENSG00EFGH': {'panels': [123]},
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
                'new': [],
                'panels': [137],
            },
            'ENSG00IJKL': {
                'symbol': 'IJKL',
                'moi': 'Mono_And_Biallelic',
                'new': [137],
                'panels': [123, 137],
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
    assert read_panels_from_participant_file(str(tmp_json)) == {1, 2, 3, 9, 99}
