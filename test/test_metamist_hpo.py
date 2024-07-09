"""
test file for metamist panel-participant matching
"""

import networkx as nx
import pytest
from obonet import read_obo
from talos.GeneratePanelData import (
    get_panels,
    match_hpo_terms,
    match_hpos_to_panels,
    match_participants_to_panels,
    update_hpo_with_label,
)
from talos.models import ParticipantHPOPanels, PhenotypeMatchedPanels

PANELAPP = 'https://fake_panelapp.agha.umccr.org/api/v1/panels/'


@pytest.fixture(name='fake_panelapp_overview')
def fixture_fake_panelapp_overview(requests_mock, fake_panelapp_overview):
    """
    a new fixture to contain the panel data
    """
    requests_mock.register_uri('GET', PANELAPP, json=fake_panelapp_overview)


def test_get_panels(fake_panelapp_overview):  # noqa: ARG001
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
    assert match_hpo_terms(panel_map=panel_map, hpo_tree=obo_parsed, hpo_str='HP:4') == {1, 2}
    assert match_hpo_terms(panel_map=panel_map, hpo_tree=obo_parsed, hpo_str='HP:2') == {1, 2}


def test_match_hpos_to_panels(fake_obo_path):
    """
    test the hpo-to-panel matching
    """
    panel_map = {'HP:2': {1, 2}, 'HP:5': {5}}
    assert match_hpos_to_panels(panel_map, fake_obo_path, all_hpos={'HP:4', 'HP:7a'}) == (
        {'HP:4': {1, 2}, 'HP:7a': {5}},
        {'HP:4': 'Goblet of Fire', 'HP:7a': 'Deathly Hallows'},
    )
    # full depth from the terminal node should capture all panels
    assert match_hpos_to_panels(panel_map, fake_obo_path, all_hpos={'HP:4', 'HP:7a'}) == (
        {'HP:4': {1, 2}, 'HP:7a': {5}},
        {'HP:4': 'Goblet of Fire', 'HP:7a': 'Deathly Hallows'},
    )


def test_read_hpo_tree(fake_obo_path):
    """
    check that reading the obo tree works
    """
    obo_parsed = read_obo(fake_obo_path)
    assert isinstance(obo_parsed, nx.MultiDiGraph)
    assert list(nx.bfs_edges(obo_parsed, 'HP:1', reverse=True)) == [
        ('HP:1', 'HP:2'),
        ('HP:2', 'HP:3'),
        ('HP:3', 'HP:4'),
        ('HP:4', 'HP:5'),
        ('HP:5', 'HP:6'),
        ('HP:6', 'HP:7a'),
        ('HP:6', 'HP:7b'),
    ]
    assert obo_parsed.nodes()['HP:3'] == {
        'name': 'Prisoner of Azkaban',
        'comment': 'where my Hippogriff at?',
        'is_a': ['HP:2'],
    }


def test_match_participants_to_panels():
    """ """
    party_hpo = PhenotypeMatchedPanels(
        samples={
            'luke_skywalker': {
                'external_id': 'participant1',
                'family_id': 'fam1',
                'hpo_terms': [{'id': 'HP:1', 'label': ''}, {'id': 'HP:2', 'label': ''}],
                'panels': {137},
            },
            'participant2': {
                'external_id': 'participant2',
                'family_id': 'fam2',
                'hpo_terms': [{'id': 'HP:1', 'label': ''}, {'id': 'HP:6', 'label': ''}],
                'panels': {137},
            },
        },
    )
    hpo_to_panels = {'HP:1': {101, 102}, 'HP:2': {2002}, 'HP:3': {00, 1, 2}, 'HP:6': {666}}
    match_participants_to_panels(participant_hpos=party_hpo, hpo_panels=hpo_to_panels)
    assert party_hpo.samples['luke_skywalker'] == ParticipantHPOPanels(
        external_id='participant1',
        family_id='fam1',
        hpo_terms=[{'id': 'HP:1', 'label': ''}, {'id': 'HP:2', 'label': ''}],
        panels={137, 101, 102, 2002},
    )
    assert party_hpo.samples['participant2'] == ParticipantHPOPanels(
        external_id='participant2',
        family_id='fam2',
        hpo_terms=[{'id': 'HP:1', 'label': ''}, {'id': 'HP:6', 'label': ''}],
        panels={137, 101, 102, 666},
    )


def test_update_hpo_with_description():
    """
    test that the description is added to the hpo
    """
    hpo_dict = PhenotypeMatchedPanels(
        samples={
            'luke_skywalker': {
                'external_id': 'participant1',
                'family_id': 'fam1',
                'hpo_terms': [{'id': 'HP:1', 'label': ''}, {'id': 'HP:2', 'label': ''}],
                'panels': {137},
            },
        },
    )
    hpo_to_desc = {'HP:1': "Philosopher's Stone", 'HP:2': 'Chamber of Secrets'}
    hpo_dict = update_hpo_with_label(hpo_dict, hpo_to_desc)
    assert hpo_dict == PhenotypeMatchedPanels(
        samples={
            'luke_skywalker': {
                'external_id': 'participant1',
                'family_id': 'fam1',
                'hpo_terms': [
                    {'id': 'HP:1', 'label': "Philosopher's Stone"},
                    {'id': 'HP:2', 'label': 'Chamber of Secrets'},
                ],
                'panels': {137},
            },
        },
    )
