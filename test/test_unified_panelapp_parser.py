from talos.UnifiedPanelAppParser import match_hpos_to_panels, match_participants_to_panels
from talos.models import PanelApp, PhenoPacketHpo, ParticipantHPOPanels, DownloadedPanelApp, PanelShort, CURRENT_VERSION


def test_match_hpos_to_panels(fake_obo_path):
    """
    test the hpo-to-panel matching
    this has now been adjusted to account for the full leaf-to-root traversal, instead of stopping at 3 layers
    """
    hpo_panel_map = {
        1: [PhenoPacketHpo(id='HP:2', label='')],
        2: [PhenoPacketHpo(id='HP:2', label='')],
        5: [PhenoPacketHpo(id='HP:5', label='')],
    }

    # full depth from the terminal node should capture all panels
    result = match_hpos_to_panels(hpo_panel_map, fake_obo_path, all_hpos={'HP:4', 'HP:7a'})

    assert result == {
        'HP:4': {1, 2},
        'HP:7a': {1, 2, 5},
    }


def test_match_participants_to_panels():
    """

    Returns:

    """

    cached_panelapp = DownloadedPanelApp(
        versions=[
            PanelShort(id=1, version='2'),
            PanelShort(id=3, version='4'),
            PanelShort(id=5, version='6'),
        ]
    )

    papp = PanelApp(
        participants={
            'sam1': ParticipantHPOPanels(
                hpo_terms=[
                    PhenoPacketHpo(id='HP:1', label='1'),
                    PhenoPacketHpo(id='HP:2', label='2'),
                ]
            ),
            'sam2': ParticipantHPOPanels(
                hpo_terms=[
                    PhenoPacketHpo(id='HP:1', label='1'),
                ]
            ),
        }
    )

    hpo_panels = {'HP:1': {1, 3}, 'HP:2': {1, 5}}

    match_participants_to_panels(
        panelapp_data=papp,
        hpo_panels=hpo_panels,
        cached_panelapp=cached_panelapp,
    )

    assert dict(papp) == {
        'metadata': [
            PanelShort(id=1, name='', version='2'),
            PanelShort(id=3, name='', version='4'),
            PanelShort(id=5, name='', version='6'),
        ],
        'genes': {},
        'participants': {
            'sam1': ParticipantHPOPanels(
                external_id='',
                family_id='',
                hpo_terms=[PhenoPacketHpo(id='HP:1', label='1'), PhenoPacketHpo(id='HP:2', label='2')],
                panels={1, 3, 5},
                matched_genes=set(),
                matched_phenotypes=set(),
            ),
            'sam2': ParticipantHPOPanels(
                external_id='',
                family_id='',
                hpo_terms=[PhenoPacketHpo(id='HP:1', label='1')],
                panels={1, 3},
                matched_genes=set(),
                matched_phenotypes=set(),
            ),
        },
        'version': CURRENT_VERSION,
    }
