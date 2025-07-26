import pytest

from talos.models import CURRENT_VERSION, DownloadedPanelApp, PanelApp, PanelShort, ParticipantHPOPanels, PhenoPacketHpo
from talos.static_values import get_granular_date
from talos.UnifiedPanelAppParser import (
    get_simple_moi,
    match_hpos_to_panels,
    match_participants_to_panels,
)


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
            # 99 is a forced panel in the test config
            PanelShort(id=99, version='99'),
        ],
    )

    papp = PanelApp(
        participants={
            'sam1': ParticipantHPOPanels(
                hpo_terms=[
                    PhenoPacketHpo(id='HP:1', label='1'),
                    PhenoPacketHpo(id='HP:2', label='2'),
                ],
            ),
            'sam2': ParticipantHPOPanels(
                hpo_terms=[
                    PhenoPacketHpo(id='HP:1', label='1'),
                ],
            ),
        },
    )

    hpo_panels = {'HP:1': {1, 3}, 'HP:2': {1, 5}}

    match_participants_to_panels(
        panelapp_data=papp,
        hpo_panels=hpo_panels,
        cached_panelapp=cached_panelapp,
    )

    assert dict(papp) == {
        'creation_date': get_granular_date(),
        'metadata': {
            1: PanelShort(id=1, name='', version='2'),
            3: PanelShort(id=3, name='', version='4'),
            5: PanelShort(id=5, name='', version='6'),
            99: PanelShort(id=99, name='', version='99'),
        },
        'genes': {},
        'participants': {
            'sam1': ParticipantHPOPanels(
                external_id='',
                family_id='',
                hpo_terms=[PhenoPacketHpo(id='HP:1', label='1'), PhenoPacketHpo(id='HP:2', label='2')],
                panels={99, 1, 3, 5},
                matched_genes=set(),
                matched_phenotypes=set(),
            ),
            'sam2': ParticipantHPOPanels(
                external_id='',
                family_id='',
                hpo_terms=[PhenoPacketHpo(id='HP:1', label='1')],
                panels={99, 1, 3},
                matched_genes=set(),
                matched_phenotypes=set(),
            ),
        },
        'version': CURRENT_VERSION,
    }


@pytest.mark.parametrize(
    'strings,expected,chrom',
    [
        (set(), 'Biallelic', '1'),
        (set(), 'Hemi_Bi_In_Female', 'X'),
        ({'blag'}, 'Biallelic', '1'),
        ({'blag'}, 'Hemi_Bi_In_Female', 'X'),
        ({'biallelic ANY'}, 'Biallelic', '1'),
        ({'both something', 'something'}, 'Mono_And_Biallelic', '1'),
        ({'monoallelic', 'something'}, 'Monoallelic', '1'),
        ({'monoallelic', 'something'}, 'Hemi_Mono_In_Female', 'X'),
        ({'monoallelic', 'biallelic'}, 'Mono_And_Biallelic', '1'),
        ({'monoallelic', 'biallelic'}, 'Hemi_Mono_In_Female', 'X'),
        ({'x-linked'}, 'Hemi_Mono_In_Female', 'X'),
        ({'x-linked biallelic'}, 'Hemi_Bi_In_Female', 'X'),
    ],
)
def test_get_simple_moi(strings: set[str], expected: str, chrom: str):
    """
    Tests the string parsing down to simple representation
    """
    assert get_simple_moi(strings, chrom) == expected
