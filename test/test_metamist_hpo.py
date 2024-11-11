"""
test file for metamist panel-participant matching
"""

from talos.GeneratePanelData import get_panels, match_hpos_to_panels, match_participants_to_panels
from talos.models import ParticipantHPOPanels, PhenotypeMatchedPanels


def test_get_panels(httpx_mock, fake_panelapp_overview):
    """
    check that the endpoint parser works ok
    """
    httpx_mock.add_response(url='https://panelapp.agha.umccr.org/api/v1/panels/', json=fake_panelapp_overview)
    panels_parsed = get_panels('https://panelapp.agha.umccr.org/api/v1/panels/')
    assert panels_parsed == {'HP:1': {2}, 'HP:4': {1}, 'HP:6': {2}}


def test_match_hpos_to_panels():
    """
    test the hpo-to-panel matching
    """
    panel_map = {'HP:0000574': {1, 2}, 'HP:0000234': {5}}
    assert match_hpos_to_panels(panel_map, all_hpos={'HP:0000574', 'HP:0000234'}) == {
        'HP:0000574': {1, 2, 5},
        'HP:0000234': {5},
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
