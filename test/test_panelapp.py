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
    combine_mendeliome_with_other_panels,
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
    with open(PANELAPP_OLDER, 'r', encoding='utf-8') as handle:
        requests_mock.register_uri(
            'GET',
            f'https://panelapp.agha.umccr.org/api/v1/panels/137?version={OLD_VERSION}',
            complete_qs=True,
            json=json.load(handle),
        )


def test_panel_query(fake_panelapp):  # pylint: disable=unused-argument
    """
    check that the default parsing delivers correct data
    :param fake_panelapp:
    :return:
    """
    result = get_panel_green(panel_id='137')
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
    main = get_panel_green('137')
    early_keys = set(main.keys())
    additional = get_panel_green('126')
    combine_mendeliome_with_other_panels(main, additional)
    assert set(main.keys()) == early_keys
    assert main['ENSG00ABCD'].get('flags') == ['Incidentalome']


def test_second_panel_update_moi():
    """
    panels overlap with MOI in second panel
    - expect overwriting of MOI
    """
    main = {
        'metadata': {'panel_name': 'NAME', 'additional_panels': []},
        'ensg1': {'entity_name': '1', 'moi': None, 'flags': []},
    }
    additional = {
        'metadata': {
            'panel_name': 'NAME',
            'panel_version': '1.0',
            'panel_id': '123',
            'additional_panels': [],
        },
        'ensg1': {'entity_name': '1', 'moi': 'REALLY_BIG'},
    }
    combine_mendeliome_with_other_panels(main, additional)
    assert main['ensg1'].get('flags') == ['NAME']
    assert main['ensg1'].get('moi') == 'REALLY_BIG'
    assert len(main['metadata']['additional_panels']) == 1
    assert main['metadata']['additional_panels'][0]['panel_id'] == '123'


def test_second_panel_no_overlap():
    """
    two panels don't overlap
    - expect both results
    """
    main = {
        'metadata': {
            'panel_name': 'NAME',
            'panel_version': '1.0',
            'panel_id': '1',
            'additional_panels': [],
        },
        'ensg1': {'entity_name': '1', 'moi': None, 'flags': []},
    }
    additional = {
        'metadata': {
            'panel_name': 'ADD',
            'panel_version': '1.1',
            'panel_id': '2',
            'additional_panels': [],
        },
        'ensg2': {'entity_name': '2', 'moi': 'REALLY_BIG'},
    }
    combine_mendeliome_with_other_panels(main, additional)
    assert main['ensg1'].get('flags') == []
    assert main['ensg1'].get('moi') is None
    assert main['ensg2'].get('flags') == ['ADD']
    assert main['ensg2'].get('moi') == 'REALLY_BIG'
    assert len(main['metadata']['additional_panels']) == 1
    assert main['metadata']['additional_panels'][0]['panel_version'] == '1.1'
