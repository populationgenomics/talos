"""
script testing methods within reanalysis/validate_categories.py
"""

from mendelbrot.pedigree_parser import PedigreeParser

from talos.models import (
    Coordinates,
    PanelApp,
    ReportVariant,
    ResultData,
    ResultMeta,
    SmallVariant,
    StructuralVariant,
)
from talos.static_values import get_granular_date
from talos.validate_moi import count_families, filter_results_to_panels, prepare_results_shell, set_up_moi_filters
from test.test_utils import ONE_EXPECTED, THREE_EXPECTED, TWO_EXPECTED, ZERO_EXPECTED

TEST_COORDS = Coordinates(chrom='1', pos=1, ref='A', alt='C')
TEST_COORDS_2 = Coordinates(chrom='2', pos=2, ref='G', alt='T')
VAR_1 = SmallVariant(coordinates=TEST_COORDS, info={}, transcript_consequences=[])
VAR_2 = SmallVariant(coordinates=TEST_COORDS_2, info={}, transcript_consequences=[])
REP_SAM1_1 = ReportVariant(sample='sam1', var_data=VAR_1, categories={'ClinvarP/LP': get_granular_date()}, gene='ENSG1')
REP_SAM3_1 = ReportVariant(sample='sam3', var_data=VAR_1, categories={'ClinvarP/LP': get_granular_date()}, gene='ENSG4')
REP_SAM3_2 = ReportVariant(sample='sam3', var_data=VAR_2, categories={'2': get_granular_date()}, gene='ENSG5')

dirty_data = [REP_SAM1_1, REP_SAM3_1, REP_SAM3_2]


def test_results_shell(pedigree_path: str):
    """

    Returns:

    """
    panelapp = PanelApp(
        metadata={1: {'id': 1, 'name': 'lorem'}, 2: {'id': 2, 'name': 'ipsum'}, 3: {'id': 3, 'name': 'etc'}},
        genes={'ENSG1': {'symbol': 'G1'}},
        participants={
            'male': {'panels': {1, 3}, 'external_id': 'male', 'hpo_terms': [{'id': 'HPB', 'label': 'Boneitis!'}]},
            'father_1': {'external_id': 'father_1', 'hpo_terms': []},
            'mother_1': {'external_id': 'mother_1', 'hpo_terms': []},
            'female': {
                'panels': {1, 2},
                'external_id': 'female',
                'hpo_terms': [{'id': 'HPF', 'label': 'HPFemale'}],
            },
            'father_2': {'external_id': 'father_2', 'hpo_terms': []},
            'mother_2': {'external_id': 'mother_2', 'hpo_terms': []},
        },
    )
    source_samples = {'small': {'male'}, 'sv': {'female'}}
    shell = prepare_results_shell(
        results_meta=ResultMeta(),
        source_samples=source_samples,
        pedigree=PedigreeParser(pedigree_path=pedigree_path),
        panelapp=panelapp,
    )

    # top level only has the two affected participants
    expected = ResultData(
        results={
            'male': {
                'metadata': {
                    'ext_id': 'male',
                    'family_id': 'family_1',
                    'members': {
                        'male': {'sex': 'male', 'affected': True},
                        'father_1': {'sex': 'male', 'affected': False},
                        'mother_1': {'sex': 'female', 'affected': False},
                    },
                    'phenotypes': [{'id': 'HPB', 'label': 'Boneitis!'}],
                    'panel_details': {1: {'id': 1, 'name': 'lorem'}, 3: {'id': 3, 'name': 'etc'}},
                    'present_in_small': True,
                },
            },
            'female': {
                'metadata': {
                    'ext_id': 'female',
                    'family_id': 'family_2',
                    'members': {
                        'female': {'sex': 'female', 'affected': True, 'ext_id': 'female'},
                        'father_2': {'sex': 'male', 'affected': False, 'ext_id': 'father_2'},
                        'mother_2': {'sex': 'female', 'affected': False, 'ext_id': 'mother_2'},
                    },
                    'phenotypes': [{'id': 'HPF', 'label': 'HPFemale'}],
                    'panel_details': {1: {'id': 1, 'name': 'lorem'}, 2: {'id': 2, 'name': 'ipsum'}},
                    'solved': True,
                    'present_in_sv': True,
                },
            },
        },
    )

    assert shell.results == expected.results


def test_gene_clean_results_no_personal(caplog):
    """
    tests the per-participant gene-filtering of results
    messy test, write and pass file paths
    """
    results_holder = ResultData(
        results={
            'sam1': {'metadata': {'ext_id': 'sam1', 'family_id': 'family_1'}},
            'sam2': {'metadata': {'ext_id': 'sam2', 'family_id': 'family_2'}},
            'sam3': {'metadata': {'ext_id': 'sam3', 'family_id': 'family_3'}},
        },
    )

    panelapp = PanelApp(
        metadata={
            137: {'id': 137, 'version': '137'},
            1: {'id': 1, 'version': '1', 'name': '1'},
            2: {'id': 2, 'version': '2', 'name': '2'},
            3: {'id': 3, 'version': '3', 'name': '3'},
            4: {'id': 4, 'version': '4', 'name': '4'},
        },
        genes={
            'ENSG1': {'panels': {137, 1}, 'symbol': 'G1'},
            'ENSG2': {'symbol': 'G2'},
            'ENSG3': {'panels': {2}, 'symbol': 'G3'},
            'ENSG4': {'panels': {3}, 'new': [3], 'symbol': 'G4'},
            'ENSG5': {'panels': {4}, 'symbol': 'G5'},
        },
    )

    filter_results_to_panels(results_holder=results_holder, result_list=dirty_data, panelapp=panelapp)

    for sample_id in ['sam1', 'sam3']:
        assert f'Participant {sample_id} not found in panelapp participants' in caplog.text

    assert len(results_holder.results['sam1'].variants) == ONE_EXPECTED
    assert results_holder.results['sam1'].variants[0].gene == 'ENSG1'
    assert results_holder.results['sam1'].variants[0].flags == set()
    assert len(results_holder.results['sam2'].variants) == ZERO_EXPECTED


def test_gene_clean_results_personal():
    """
    tests the per-participant gene-filtering of results
    messy test, write and pass file paths
    """
    results_holder = ResultData(
        results={
            'sam1': {'metadata': {'ext_id': 'sam1', 'family_id': 'family_1', 'panel_ids': {1}}},
            'sam2': {'metadata': {'ext_id': 'sam2', 'family_id': 'family_2'}},
            'sam3': {'metadata': {'ext_id': 'sam3', 'family_id': 'family_3', 'panel_ids': {1}}},
        },
    )
    panelapp = PanelApp(
        metadata={
            137: {'id': 137, 'version': '137'},
            1: {'id': 1, 'version': '1', 'name': '1'},
            2: {'id': 2, 'version': '2', 'name': '2'},
            3: {'id': 3, 'version': '3', 'name': '3'},
            4: {'id': 4, 'version': '4', 'name': '4'},
        },
        genes={
            'ENSG1': {'panels': {137, 1}, 'symbol': 'G1'},
            'ENSG2': {'symbol': 'G2'},
            'ENSG3': {'panels': {2}, 'symbol': 'G3'},
            'ENSG4': {'panels': {3}, 'new': [3], 'symbol': 'G4'},
            'ENSG5': {'panels': {4}, 'symbol': 'G5'},
        },
        participants={
            'sam1': {'panels': {1}, 'hpo_terms': [{'id': 'HP1', 'label': 'HP1'}]},
            'sam2': {'hpo_terms': [{'id': 'HP2', 'label': 'HP2'}]},
            'sam3': {'panels': {3, 4}, 'hpo_terms': [{'id': 'HP3', 'label': 'HP3'}]},
        },
    )

    filter_results_to_panels(results_holder, dirty_data, panelapp=panelapp)
    assert len(results_holder.results['sam1'].variants) == ONE_EXPECTED
    assert results_holder.results['sam1'].variants[0].gene == 'ENSG1'
    assert not results_holder.results['sam1'].variants[0].flags
    assert results_holder.results['sam1'].variants[0].panels.matched == {1: '1'}
    assert not results_holder.results['sam2'].variants
    assert len(results_holder.results['sam3'].variants) == TWO_EXPECTED
    for event in results_holder.results['sam3'].variants:
        if event.gene == 'ENSG4':
            assert event.panels.matched == {3: '3'}
        if event.gene == 'ENSG5':
            assert event.panels.matched == {4: '4'}


def test_gene_clean_results_max_confidence():
    """
    max_confidence is the highest confidence across every panel applied to the participant which carries the gene
    - the custom panel (0) is never assigned to a participant, it must arrive through the forced panel set
    - a panel with no recorded confidence contributes 0
    - a participant missing from panelapp gives 0, not a crash
    """
    results_holder = ResultData(
        results={
            'sam1': {'metadata': {'ext_id': 'sam1', 'family_id': 'family_1'}},
            'sam2': {'metadata': {'ext_id': 'sam2', 'family_id': 'family_2'}},
        },
    )
    panelapp = PanelApp(
        metadata={
            0: {'id': 0, 'name': 'custom'},
            137: {'id': 137, 'version': '137'},
            1: {'id': 1, 'version': '1', 'name': '1'},
            2: {'id': 2, 'version': '2', 'name': '2'},
        },
        genes={
            # on the custom panel and one participant panel - custom wins
            'ENSG_CUSTOM': {'panels': {0, 1}, 'panel_confidences': {0: 3, 1: 2}, 'symbol': 'GC'},
            # only on participant panels, take the highest of those
            'ENSG_AMBER': {'panels': {137, 1}, 'panel_confidences': {137: 1, 1: 2}, 'symbol': 'GA'},
            # on a participant panel, but no confidence recorded for it
            'ENSG_NONE': {'panels': {2}, 'symbol': 'GN'},
        },
        participants={
            'sam1': {'panels': {137, 1, 2}, 'hpo_terms': []},
        },
    )
    events = [
        ReportVariant(sample='sam1', var_data=VAR_1, categories={'2': get_granular_date()}, gene='ENSG_CUSTOM'),
        ReportVariant(sample='sam1', var_data=VAR_2, categories={'2': get_granular_date()}, gene='ENSG_AMBER'),
        ReportVariant(sample='sam1', var_data=VAR_1, categories={'2': get_granular_date()}, gene='ENSG_NONE'),
        # sam2 is not in panelapp.participants, but the gene is on the custom panel so it survives
        ReportVariant(sample='sam2', var_data=VAR_1, categories={'2': get_granular_date()}, gene='ENSG_CUSTOM'),
    ]

    filter_results_to_panels(results_holder, events, panelapp=panelapp)

    by_gene = {event.gene: event for event in results_holder.results['sam1'].variants}
    assert by_gene['ENSG_CUSTOM'].max_confidence == THREE_EXPECTED
    assert by_gene['ENSG_CUSTOM'].panels.forced == {0: 'custom'}
    assert by_gene['ENSG_AMBER'].max_confidence == TWO_EXPECTED
    assert by_gene['ENSG_NONE'].max_confidence == ZERO_EXPECTED

    # no participant panels at all, but the custom panel is forced, so it still counts as green
    assert len(results_holder.results['sam2'].variants) == ONE_EXPECTED
    assert results_holder.results['sam2'].variants[0].max_confidence == THREE_EXPECTED


def test_update_results_meta(pedigree_path: str):
    """
    testing the dict update
    """
    pedigree = PedigreeParser(pedigree_path)
    ped_samples = {'male', 'female', 'mother_1', 'father_1', 'mother_2', 'father_2'}
    pedigree.set_participants(pedigree.strip_pedigree_to_samples(ped_samples))

    assert count_families(pedigree=pedigree) == {
        'affected': TWO_EXPECTED,
        'male': THREE_EXPECTED,
        'female': THREE_EXPECTED,
        'trios': TWO_EXPECTED,
    }


def test_count_families_missing_father(pedigree_path: str):
    """
    testing the dict update
    """

    assert count_families(pedigree=PedigreeParser(pedigree_path)) == {
        'affected': TWO_EXPECTED,
        'male': THREE_EXPECTED,
        'female': THREE_EXPECTED,
        'trios': TWO_EXPECTED,
    }


def test_count_families_quad(quad_ped: str):
    """
    testing the dict update
    this one is being marked as a trio as there are 4 samples, but only one affected
    my vibe here is that a 'quad' is typically 2 children 2 parents
    for affected child, unaffected sib, 2 parents... that's just a trio + 1?
    """

    ped_samples = {'PROBAND', 'SIBLING', 'FATHER', 'MOTHER'}
    ped = PedigreeParser(quad_ped)
    ped.set_participants(ped.strip_pedigree_to_samples(ped_samples))
    assert count_families(pedigree=ped) == {'affected': 1, 'male': 3, 'female': 1, 'trios': 1}


def test_regions_in_moi_filters_and_panel_filtering(pedigree_path: str):
    """region MOIs get a MOI runner, and region events are matched to panels via the region entry"""
    panelapp = PanelApp(
        metadata={137: {'id': 137, 'name': 'Mendeliome'}},
        genes={'ENSG1': {'symbol': 'G1', 'moi': 'Biallelic'}},
        regions={'REGION': {'name': 'REGION', 'moi': 'Monoallelic', 'panels': {137}, 'panel_confidences': {137: 3}}},
        participants={'sam1': {'panels': {137}}},
    )
    assert set(set_up_moi_filters(panelapp, PedigreeParser(pedigree_path=pedigree_path))) == {
        'Biallelic',
        'Monoallelic',
    }

    sv = StructuralVariant(coordinates=TEST_COORDS, info={'gene_id': 'REGION'})
    results_holder = ResultData(results={'sam1': {'metadata': {'family_id': 'family_1'}}})
    filter_results_to_panels(
        results_holder=results_holder,
        result_list=[ReportVariant(sample='sam1', var_data=sv, gene='REGION')],
        panelapp=panelapp,
    )
    assert len(results_holder.results['sam1'].variants) == ONE_EXPECTED
    assert results_holder.results['sam1'].variants[0].max_confidence == 3
