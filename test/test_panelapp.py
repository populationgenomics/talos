"""
tests for the PanelApp parser
"""

import json

import pytest

from reanalysis.query_panelapp import (
    get_best_moi,
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


def test_panel_query(fake_panelapp):  # pylint: disable=unused-argument
    """
    check that the default parsing delivers correct data
    :param fake_panelapp: fake web hook mock
    """

    gd = {'genes': {}, 'metadata': []}
    old_data = {'ENSG00ABCD': [1], 'ENSG00EFGH': [137]}
    get_panel_green(gd, old_data=old_data)
    assert gd['genes']['ENSG00ABCD']['moi'] == {'biallelic'}
    assert gd['genes']['ENSG00ABCD']['panels'] == [137]
    assert gd['genes']['ENSG00EFGH']['moi'] == {'monoallelic'}
    assert old_data['ENSG00ABCD'] == [1, 137]


def test_panel_query_addition(fake_panelapp):  # pylint: disable=unused-argument
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
                'moi': {'monoallelic'},
                'new': [],
                'panels': [137],
            },
            'ENSG00IJKL': {
                'symbol': 'IJKL',
                'moi': {'both'},
                'new': [137],
                'panels': [123, 137],
            },
        },
    }

    # should query for and integrate the incidentalome content
    get_panel_green(
        gd, panel_id=126, old_data={'ENSG00EFGH': [137, 126], 'ENSG00IJKL': [137]}
    )
    assert gd['genes']['ENSG00ABCD']['moi'] == {'monoallelic', 'biallelic'}
    assert gd['genes']['ENSG00ABCD']['panels'] == [137, 126]
    assert gd['genes']['ENSG00IJKL']['moi'] == {'both'}
    assert gd['genes']['ENSG00IJKL']['panels'] == [123, 137]
    assert 'ENSG00EFGH' not in gd['genes']


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


def test_get_best_moi_empty():
    """
    check that the MOI summary works
    """

    d = {'ensg1': {'moi': {}, 'chrom': '1'}}
    get_best_moi(d)
    assert d['ensg1']['moi'] == 'Biallelic'

    d = {'ensg1': {'moi': {}, 'chrom': 'X'}}
    get_best_moi(d)
    assert d['ensg1']['moi'] == 'Hemi_Bi_In_Female'


def test_get_best_moi_mono():
    """
    check that the MOI summary works
    """

    d = {'ensg1': {'moi': {'monoallelic'}, 'chrom': '1'}}
    get_best_moi(d)
    assert d['ensg1']['moi'] == 'Monoallelic'


def test_get_best_moi_mono_and_biallelic():
    """
    check that the MOI summary works
    """

    d = {'ensg1': {'moi': {'monoallelic', 'biallelic'}, 'chrom': '1'}}
    get_best_moi(d)
    assert d['ensg1']['moi'] == 'Mono_And_Biallelic'


def test_get_best_moi_1():
    """
    check that the MOI summary works
    """

    d = {'ensg1': {'moi': {'Monoallelic', 'Biallelic', 'both'}, 'chrom': '1'}}
    get_best_moi(d)
    assert d['ensg1']['moi'] == 'Mono_And_Biallelic'


def test_get_best_moi_x():
    """
    check that the MOI summary works
    """

    d = {'ensg1': {'moi': {'x-linked biallelic', 'x-linked'}, 'chrom': 'X'}}
    get_best_moi(d)
    assert d['ensg1']['moi'] == 'Hemi_Mono_In_Female'
