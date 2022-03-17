"""
tests for the PanelApp parser
"""

import json
import os

import pytest

from reanalysis.query_panelapp import (
    get_json_response,
    get_panel_green,
    get_panel_changes,
)

PWD = os.path.dirname(__file__)
INPUT = os.path.join(PWD, 'input')
PANELAPP_LATEST = os.path.join(INPUT, 'panelapp_current_137.json')
PANELAPP_OLDER = os.path.join(INPUT, 'panelapp_older_137.json')
LATEST_EXPECTED = os.path.join(INPUT, 'panel_green_latest_expected.json')
CHANGES_EXPECTED = os.path.join(INPUT, 'panel_changes_expected.json')
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
            'https://panelapp.agha.umccr.org/api/v1/137',
            json=json.load(handle),
        )
    with open(PANELAPP_OLDER, 'r', encoding='utf-8') as handle:
        requests_mock.register_uri(
            'GET',
            f'https://panelapp.agha.umccr.org/api/v1/137?version={OLD_VERSION}',
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
        assert result == json.load(handle)


def test_get_panel_changes(fake_panelapp):  # pylint: disable=unused-argument
    """

    :param fake_panelapp:
    :return:
    """
    with open(LATEST_EXPECTED, 'r', encoding='utf-8') as handle:
        latest = json.load(handle)

    get_panel_changes(
        previous_version=OLD_VERSION, panel_id='137', latest_content=latest
    )
    with open(CHANGES_EXPECTED, 'r', encoding='utf-8') as handle2:
        assert latest == json.load(handle2)


def test_get_json_response(fake_panelapp):  # pylint: disable=unused-argument
    """
    read the json content via an API call
    :param fake_panelapp:
    :return:
    """
    result = get_json_response('https://panelapp.agha.umccr.org/api/v1/137')
    with open(PANELAPP_LATEST, 'r', encoding='utf-8') as handle:
        assert result == json.load(handle)
