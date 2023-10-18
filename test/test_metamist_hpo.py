"""
test file for metamist panel-participant matching
"""

import pytest

from obonet import read_obo

from helpers.hpo_panel_match import (
    match_hpos_to_panels,
    match_participants_to_panels,
    get_panels,
    match_hpo_terms,
)

PANELAPP = 'https://fake_panelapp.agha.umccr.org/api/v1/panels/'


@pytest.fixture(name='fake_panelapp_overview')
def fixture_fake_panelapp_overview(requests_mock, fake_panelapp_overview):
    """
    a new fixture to contain the panel data
    :param requests_mock:
    :param fake_panelapp_overview:
    """
    requests_mock.register_uri(
        'GET',
        PANELAPP,
        json=fake_panelapp_overview,
    )


def test_get_panels(fake_panelapp_overview):  # pylint: disable=unused-argument
    """
    check that the endpoint parser works ok
    """
    panels_parsed = get_panels(PANELAPP)
    assert panels_parsed == {'HP:1': {2}, 'HP:4': {1}, 'HP:6': {2}}


def test_match_hpo_terms(fake_obo_path):
    """
    check that HP tree traversal works
    this test is kinda limited now that the layer count is constant
    """
    obo_parsed = read_obo(fake_obo_path)
    panel_map = {'HP:2': {1, 2}}
    assert match_hpo_terms(
        panel_map=panel_map, hpo_tree=obo_parsed, hpo_str='HP:4'
    ) == {1, 2}
    assert match_hpo_terms(
        panel_map=panel_map, hpo_tree=obo_parsed, hpo_str='HP:2'
    ) == {1, 2}
    # assert (
    #     match_hpo_terms(panel_map=panel_map, hpo_tree=obo_parsed, hpo_str='HP:3')
    #     == set()
    # )


def test_match_hpos_to_panels(fake_obo_path):
    """
    test the hpo-to-panel matching
    """
    panel_map = {'HP:2': {1, 2}, 'HP:5': {5}}
    assert match_hpos_to_panels(
        panel_map, fake_obo_path, all_hpos={'HP:4', 'HP:7a'}
    ) == {
        'HP:4': {1, 2},
        'HP:7a': {5},
    }
    # full depth from the terminal node should capture all panels
    assert match_hpos_to_panels(
        panel_map, fake_obo_path, all_hpos={'HP:4', 'HP:7a'}
    ) == {
        'HP:4': {1, 2},
        'HP:7a': {5},
    }


def test_match_participants_to_panels():
    """

    Returns
    -------

    """
    party_hpo = {
        'luke_skywalker': {
            'external_id': 'participant1',
            'family_id': 'fam1',
            'hpo_terms': {'HP:1', 'HP:2'},
            'panels': {'137'},
        },
        'participant2': {
            'external_id': 'participant2',
            'family_id': 'fam2',
            'hpo_terms': {'HP:1', 'HP:6'},
            'panels': {'137'},
        },
    }
    hpo_to_panels = {
        'HP:1': {'room', '101'},
        'HP:2': {'2002'},
        'HP:3': {'nothing', 'at', 'all'},
        'HP:6': {'666'},
    }
    match_participants_to_panels(participant_hpos=party_hpo, hpo_panels=hpo_to_panels)
    assert party_hpo == {
        'luke_skywalker': {
            'external_id': 'participant1',
            'family_id': 'fam1',
            'hpo_terms': sorted(['HP:1', 'HP:2']),
            'panels': sorted(['137', 'room', '101', '2002']),
        },
        'participant2': {
            'external_id': 'participant2',
            'family_id': 'fam2',
            'hpo_terms': sorted(['HP:1', 'HP:6']),
            'panels': sorted(['137', 'room', '101', '666']),
        },
    }
