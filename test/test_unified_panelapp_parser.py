from unittest.mock import MagicMock

import pytest
from obonet import read_obo

from talos.models import (
    CURRENT_VERSION,
    DownloadedPanelApp,
    HpoTerm,
    PanelApp,
    PanelDetail,
    PanelShort,
    ParticipantHPOPanels,
)
from talos.static_values import get_granular_date
from talos.UnifiedPanelAppParser import (
    CUSTOM_PANEL_ID,
    ORDERED_MOIS,
    extract_participant_data_from_pedigree,
    get_simple_moi,
    match_hpos_to_panels,
    match_participants_to_panels,
    remove_blacklisted_genes,
    update_moi_from_config,
)


def test_match_hpos_to_panels(fake_obo_path):
    """
    test the hpo-to-panel matching
    this has now been adjusted to account for the full leaf-to-root traversal, instead of stopping at 3 layers
    """
    hpo_panel_map = {
        1: [HpoTerm(id='HP:2', label='')],
        2: [HpoTerm(id='HP:2', label='')],
        5: [HpoTerm(id='HP:5', label='')],
    }

    hpo_graph = read_obo(fake_obo_path)
    # full depth from the terminal node should capture all panels
    result = match_hpos_to_panels(hpo_panel_map, hpo_graph=hpo_graph, all_hpos={'HP:4', 'HP:7a'})

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
                    HpoTerm(id='HP:1', label='1'),
                    HpoTerm(id='HP:2', label='2'),
                ],
            ),
            'sam2': ParticipantHPOPanels(
                hpo_terms=[
                    HpoTerm(id='HP:1', label='1'),
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
                hpo_terms=[HpoTerm(id='HP:1', label='1'), HpoTerm(id='HP:2', label='2')],
                panels={99, 1, 3, 5},
                matched_genes=set(),
                matched_phenotypes=set(),
            ),
            'sam2': ParticipantHPOPanels(
                external_id='',
                family_id='',
                hpo_terms=[HpoTerm(id='HP:1', label='1')],
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


def test_extract_participant_data_from_pedigree():
    pedigree = MagicMock()
    pedigree.participants = {
        'S1': MagicMock(sample_id='S1', family_id='F1', hpo_terms={'HP:1', 'HP:2'}),
        'S2': MagicMock(sample_id='S2', family_id='F1', hpo_terms={'HP:2'}),
    }
    hpo_lookup = {'HP:1': 'foo', 'HP:2': 'bar'}
    shell, all_hpos = extract_participant_data_from_pedigree(pedigree, hpo_lookup)
    assert isinstance(shell, PanelApp)
    assert all_hpos == {'HP:1', 'HP:2'}
    assert shell.participants['S1'].hpo_terms[0].label in {'foo', 'bar'}


def test_match_hpos_to_panels_no_graph():
    hpo_panel_map = {1: [HpoTerm(id='HP:1', label='foo')]}
    all_hpos = {'HP:1'}
    result = match_hpos_to_panels(hpo_panel_map, all_hpos, None)
    assert result == {}


def test_update_moi_from_config_add_new_gene():
    panelapp_data = PanelApp()
    add_genes = [{'ensg': 'ENSG1', 'symbol': 'GENE1', 'moi': ORDERED_MOIS[0], 'chrom': '1'}]
    update_moi_from_config(panelapp_data, add_genes)
    assert 'ENSG1' in panelapp_data.genes
    assert panelapp_data.genes['ENSG1'].moi == ORDERED_MOIS[0]
    assert CUSTOM_PANEL_ID in panelapp_data.genes['ENSG1'].panels


def test_update_moi_from_config_update_existing_gene():
    panelapp_data = PanelApp()
    panelapp_data.genes['ENSG2'] = PanelDetail(symbol='GENE2', moi=ORDERED_MOIS[1], panels={1}, chrom='2')
    add_genes = [{'ensg': 'ENSG2', 'moi': ORDERED_MOIS[2]}]
    update_moi_from_config(panelapp_data, add_genes)
    assert panelapp_data.genes['ENSG2'].moi == ORDERED_MOIS[2]
    assert CUSTOM_PANEL_ID in panelapp_data.genes['ENSG2'].panels


def test_remove_blacklisted_genes():
    panelapp_data = PanelApp()
    panelapp_data.genes = {'ENSG1': PanelDetail(symbol='GENE1', moi='Monoallelic', panels={1}, chrom='1')}
    remove_blacklisted_genes(panelapp_data, {'ENSG1'})
    assert 'ENSG1' not in panelapp_data.genes
