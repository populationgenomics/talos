"""
tests for the PanelApp parser
"""

from datetime import datetime
import json
import os

import pytest

from reanalysis.query_panelapp import (
    ACTIVITIES,
    get_json_response,
    get_previous_version,
    get_panel_green,
    get_panel_changes,
    PANEL_CONTENT,
)

PWD = os.path.dirname(__file__)
INPUT = os.path.join(PWD, 'input')
PANELAPP_LATEST = os.path.join(INPUT, 'panelapp_current_137.json')
PANELAPP_OLDER = os.path.join(INPUT, 'panelapp_older_137.json')
PANELAPP_ACTIVITIES = os.path.join(INPUT, 'panel_activities.json')
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
            PANEL_CONTENT.format(panel_id='137'),
            json=json.load(handle),
        )
    with open(PANELAPP_ACTIVITIES, 'r', encoding='utf-8') as handle:
        requests_mock.register_uri(
            'GET',
            ACTIVITIES.format(panel_id='137'),
            json=json.load(handle),
        )
    with open(PANELAPP_OLDER, 'r', encoding='utf-8') as handle:
        requests_mock.register_uri(
            'GET',
            f'{PANEL_CONTENT.format(panel_id="137")}?version={OLD_VERSION}',
            complete_qs=True,
            json=json.load(handle),
        )


@pytest.mark.parametrize(
    'date,version',
    [
        (datetime(year=2020, month=1, day=1), '0.2'),
        (datetime(year=2020, month=12, day=1), '0.3'),
        (datetime(year=2021, month=1, day=1), '0.3'),
        (datetime(year=1970, month=1, day=1), '0.0'),
    ],
)
def test_panel_activities(
    fake_panelapp, date, version
):  # pylint: disable=unused-argument
    """

    :param fake_panelapp:
    :return:
    """
    result = get_previous_version(panel_id='137', since=date)
    assert result == version


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
    result = get_json_response(PANEL_CONTENT.format(panel_id='137'))
    with open(PANELAPP_LATEST, 'r', encoding='utf-8') as handle:
        assert result == json.load(handle)
