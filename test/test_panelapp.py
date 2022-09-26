"""
tests for the PanelApp parser
"""

import json
import os
import pytest

from reanalysis.query_panelapp import (
    gene_list_differences,
    get_json_response,
    get_panel_green,
    parse_gene_list,
    combine_panels,
    grab_genes_only,
    get_list_from_participants,
    INCIDENTALOME,
)


PWD = os.path.dirname(__file__)
INPUT = os.path.join(PWD, 'input')
PANELAPP_LATEST = os.path.join(INPUT, 'panelapp_current_137.json')
PANELAPP_INCIDENTALOME = os.path.join(INPUT, 'incidentalome.json')
LATEST_EXPECTED = os.path.join(INPUT, 'panel_green_latest_expected.json')
CHANGES_EXPECTED = os.path.join(INPUT, 'panel_changes_expected.json')
FAKE_GENE_LIST = os.path.join(INPUT, 'fake_gene_list.json')


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


def test_panel_query(fake_panelapp):  # pylint: disable=unused-argument
    """
    check that the default parsing delivers correct data
    :param fake_panelapp:
    :return:
    """
    result = get_panel_green()
    with open(LATEST_EXPECTED, 'r', encoding='utf-8') as handle:
        expected = json.load(handle)

    assert result == expected


def test_gene_list_changes():
    """
    check usage of a gene list
    :return:
    """

    with open(LATEST_EXPECTED, 'r', encoding='utf-8') as handle:
        latest = json.load(handle)
    previous_genes = {'ABCD'}
    gene_list_differences(latest, previous_genes)
    assert not latest['ENSG00ABCD']['new']
    assert latest['ENSG00IJKL']['new']


def test_get_json_response(fake_panelapp):  # pylint: disable=unused-argument
    """
    read the json content via an API call
    :param fake_panelapp:
    :return:
    """
    result = get_json_response('https://panelapp.agha.umccr.org/api/v1/panels/137')
    with open(PANELAPP_LATEST, 'r', encoding='utf-8') as handle:
        assert result == json.load(handle)


def test_tagged_incidentalome(fake_panelapp):  # pylint: disable=unused-argument
    """
    only grab the cardiac flagged variants
    :param fake_panelapp:
    :return:
    """
    results = get_panel_green(panel_id=INCIDENTALOME, keep_tags=['cardiac'])
    key_set = set(list(results.keys()))
    assert len(key_set) == 2
    assert key_set == {'metadata', 'ENSG00WXYZ'}


def test_parse_local_gene_list():
    """
    test parsing of a local file
    :return:
    """
    found_genes = parse_gene_list(FAKE_GENE_LIST)
    assert found_genes == {'foo bar', 'foo', 'bar'}


def test_second_panel_update(fake_panelapp):  # pylint: disable=unused-argument
    """
    check that the combination method works
    """
    main = get_panel_green()
    early_keys = set(main.keys())
    additional = get_panel_green('126')
    combine_panels(main, additional)
    assert set(main.keys()) == early_keys.union({'ENSG00WXYZ'})
    assert main['ENSG00ABCD'].get('flags') == ['Incidentalome']


def test_second_panel_update_moi():
    """
    panels overlap with MOI in second panel
    - expect overwriting of MOI
    """
    main = {
        'metadata': [{'name': 'NAME', 'id': 'fish'}],
        'ensg1': {'entity_name': '1', 'moi': None, 'flags': []},
    }
    additional = {
        'metadata': [{'name': 'NAME', 'version': '1.0', 'id': '123'}],
        'ensg1': {'entity_name': '1', 'moi': 'REALLY_BIG'},
    }
    combine_panels(main, additional)
    assert main['ensg1'].get('flags') == ['NAME']
    assert main['ensg1'].get('moi') == 'REALLY_BIG'
    assert len(main['metadata']) == 2
    assert main['metadata'][0]['id'] == 'fish'
    assert main['metadata'][1]['id'] == '123'


def test_second_panel_no_overlap():
    """
    two panels don't overlap
    - expect both results
    """
    main = {
        'metadata': [{'name': 'NAME', 'version': '1.0', 'id': '1'}],
        'ensg1': {'entity_name': '1', 'moi': None, 'flags': []},
    }
    additional = {
        'metadata': [{'name': 'ADD', 'version': '1.1', 'id': '2'}],
        'ensg2': {'entity_name': '2', 'moi': 'REALLY_BIG'},
    }
    combine_panels(main, additional)
    assert main['ensg1'].get('flags') == []
    assert main['ensg1'].get('moi') is None
    assert main['ensg2'].get('flags') == ['ADD']
    assert main['ensg2'].get('moi') == 'REALLY_BIG'
    assert len(main['metadata']) == 2
    assert main['metadata'][0]['version'] == '1.0'
    assert main['metadata'][1]['version'] == '1.1'


def test_genes_only():
    """
    tests the method for getting genes from panel data
    """
    panel_data = {
        'metadata': [{'name': 'ADD', 'version': '1.1', 'id': '2'}],
        'ensg2': {'entity_name': '2', 'moi': 'REALLY_BIG'},
        'ensg3': {'entity_name': '3', 'moi': 'REALLY_LITTLE'},
        'ensgn': {'entity_name': 'n', 'moi': 'CARDBOARD_BOX'},
    }
    assert grab_genes_only(panel_data) == ['ensg2', 'ensg3', 'ensgn']


def test_get_list_from_participants():
    """
    tests the unique panel finder
    """
    party_data = {
        'i': {'panels': [1, 2], 'what': 'does'},
        'am': {'panels': [1, 3], 'the': 'fox'},
        'sam': {'panels': [9, 99], 'say?': 'Wa-pa-pa-pa-pa-pa-pow!'},
    }
    assert get_list_from_participants(party_data) == {1, 2, 3, 9, 99}
