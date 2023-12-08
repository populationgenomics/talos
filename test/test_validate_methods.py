"""
script testing methods within reanalysis/validate_categories.py
"""


from reanalysis.models import (
    Coordinates,
    PanelApp,
    PhenotypeMatchedPanels,
    ResultData,
    ResultMeta,
    # ReportPanel,
    ReportVariant,
    SmallVariant,
)
from reanalysis.validate_categories import (
    clean_and_filter,
    count_families,
    prepare_results_shell,
)


TEST_COORDS = Coordinates(chrom='1', pos=1, ref='A', alt='C')
TEST_COORDS_2 = Coordinates(chrom='2', pos=2, ref='G', alt='T')
VAR_1 = SmallVariant(coordinates=TEST_COORDS, info={}, transcript_consequences=[])
VAR_2 = SmallVariant(coordinates=TEST_COORDS_2, info={}, transcript_consequences=[])
REP_SAM1_1 = ReportVariant(
    sample='sam1', var_data=VAR_1, categories={'1'}, gene='ENSG1'
)
REP_SAM3_1 = ReportVariant(
    sample='sam3', var_data=VAR_1, categories={'1'}, gene='ENSG4'
)
REP_SAM3_2 = ReportVariant(
    sample='sam3', var_data=VAR_2, categories={'2'}, gene='ENSG5'
)


dirty_data = [REP_SAM1_1, REP_SAM3_1, REP_SAM3_2]
panel_genes = PanelApp(
    **{
        'metadata': [
            {'id': 137, 'version': '137'},
            {'id': 1, 'version': '1', 'name': '1'},
            {'id': 2, 'version': '2', 'name': '2'},
            {'id': 3, 'version': '3', 'name': '3'},
            {'id': 4, 'version': '4', 'name': '4'},
        ],
        'genes': {
            'ENSG1': {'panels': {137, 1}, 'symbol': 'G1'},
            'ENSG2': {'symbol': 'G2'},
            'ENSG3': {'panels': {2}, 'symbol': 'G3'},
            'ENSG4': {'panels': {3}, 'new': [3], 'symbol': 'G4'},
            'ENSG5': {'panels': {4}, 'symbol': 'G5'},
        },
    }
)


def test_results_shell(peddy_ped):
    """

    Returns:

    """
    samples = ['male', 'female', 'irrelevant']
    sample_panels = PhenotypeMatchedPanels(
        **{
            'samples': {
                'male': {
                    'panels': {1, 3},
                    'external_id': 'MALE!',
                    'hpo_terms': {'Boneitis'},
                },
                'female': {
                    'panels': {1, 2},
                    'external_id': 'FEMALE!',
                    'hpo_terms': {'HPfemale'},
                },
            },
            'all_panels': {1, 2, 3},
        }
    )
    panelapp = PanelApp(
        **{
            'metadata': [
                {'id': 1, 'name': 'lorem'},
                {'id': 2, 'name': 'ipsum'},
                {'id': 3, 'name': 'etc'},
            ],
            'genes': {'ENSG1': {'symbol': 'G1'}},
        }
    )
    result_meta = ResultMeta(input_file='', cohort='cohort')
    shell = prepare_results_shell(
        results_meta=result_meta,
        vcf_samples=samples,
        pedigree=peddy_ped,
        panel_data=sample_panels,
        panelapp=panelapp,
        dataset='cohort',
    )

    # top level only has the two affected participants
    expected = ResultData(
        **{
            'results': {
                'male': {
                    'metadata': {
                        'ext_id': 'MALE!',
                        'family_id': 'family_1',
                        'members': {
                            'male': {
                                'sex': 'male',
                                'affected': True,
                                'ext_id': 'MALE!',
                            },
                            'father_1': {
                                'sex': 'male',
                                'affected': False,
                                'ext_id': 'father_1',
                            },
                            'mother_1': {
                                'sex': 'female',
                                'affected': False,
                                'ext_id': 'mother_1',
                            },
                        },
                        'phenotypes': ['Boneitis'],
                        'panel_ids': [1, 3],
                        'panel_names': ['lorem', 'etc'],
                    },
                },
                'female': {
                    'metadata': {
                        'ext_id': 'FEMALE!',
                        'family_id': 'family_2',
                        'members': {
                            'female': {
                                'sex': 'female',
                                'affected': True,
                                'ext_id': 'FEMALE!',
                            },
                            'father_2': {
                                'sex': 'male',
                                'affected': False,
                                'ext_id': 'father_2',
                            },
                            'mother_2': {
                                'sex': 'female',
                                'affected': False,
                                'ext_id': 'mother_2',
                            },
                        },
                        'phenotypes': ['HPfemale'],
                        'panel_ids': [1, 2],
                        'panel_names': ['lorem', 'ipsum'],
                        'solved': True,
                    },
                },
            }
        }
    )

    assert shell.results == expected.results


def test_gene_clean_results_no_personal():
    """
    tests the per-participant gene-filtering of results
    messy test, write and pass file paths
    """
    results_holder = ResultData(
        **{
            'results': {
                'sam1': {'metadata': {'ext_id': 'sam1', 'family_id': 'family_1'}},
                'sam2': {'metadata': {'ext_id': 'sam2', 'family_id': 'family_2'}},
                'sam3': {'metadata': {'ext_id': 'sam3', 'family_id': 'family_3'}},
            }
        }
    )

    clean = clean_and_filter(
        results_holder=results_holder,
        result_list=dirty_data,
        panelapp_data=panel_genes,
        dataset='cohort',
    )
    assert len(clean.results['sam1'].variants) == 1
    assert clean.results['sam1'].variants[0].gene == 'ENSG1'
    assert clean.results['sam1'].variants[0].flags == set()
    assert len(clean.results['sam2'].variants) == 0
    assert len(clean.results['sam3'].variants) == 2
    assert {x.gene for x in clean.results['sam3'].variants} == {'ENSG4', 'ENSG5'}


def test_gene_clean_results_personal():
    """
    tests the per-participant gene-filtering of results
    messy test, write and pass file paths
    """
    results_holder = ResultData(
        **{
            'results': {
                'sam1': {
                    'metadata': {
                        'ext_id': 'sam1',
                        'family_id': 'family_1',
                        'panel_ids': {1},
                    },
                },
                'sam2': {
                    'metadata': {'ext_id': 'sam2', 'family_id': 'family_2'},
                },
                'sam3': {
                    'metadata': {
                        'ext_id': 'sam3',
                        'family_id': 'family_3',
                        'panel_ids': {1},
                    },
                },
            }
        }
    )
    personal_panels = PhenotypeMatchedPanels(
        **{
            'samples': {
                'sam1': {'panels': {1}, 'hpo_terms': {'HP1'}},
                'sam2': {'hpo_terms': {'HP2'}},
                'sam3': {'panels': {3, 4}, 'hpo_terms': {'HP3'}},
            }
        }
    )

    clean = clean_and_filter(
        results_holder=results_holder,
        result_list=dirty_data,
        panelapp_data=panel_genes,
        participant_panels=personal_panels,
        dataset='cohort',
    )
    assert len(clean.results['sam1'].variants) == 1
    assert clean.results['sam1'].variants[0].gene == 'ENSG1'
    assert not clean.results['sam1'].variants[0].flags
    assert clean.results['sam1'].variants[0].panels.matched == {'1'}
    assert not clean.results['sam2'].variants
    assert len(clean.results['sam3'].variants) == 2
    for event in clean.results['sam3'].variants:
        if event.gene == 'ENSG4':
            assert event.panels.matched == {'3'}
        if event.gene == 'ENSG5':
            assert event.panels.matched == {'4'}


def test_update_results_meta(peddy_ped):
    """
    testing the dict update
    """

    ped_samples = ['male', 'female', 'mother_1', 'father_1', 'mother_2', 'father_2']

    assert count_families(pedigree=peddy_ped, samples=ped_samples) == {
        'affected': 2,
        'male': 3,
        'female': 3,
        'trios': 2,
    }


def test_count_families_missing_father(peddy_ped):
    """
    testing the dict update
    """

    ped_samples = ['male', 'female', 'mother_1', 'mother_2', 'father_2']

    assert count_families(pedigree=peddy_ped, samples=ped_samples) == {
        'affected': 2,
        'male': 2,
        'female': 3,
        'trios': 1,
        '2': 1,
    }


def test_count_families_quad(quad_ped):
    """
    testing the dict update
    """

    ped_samples = ['PROBAND', 'SIBLING', 'FATHER', 'MOTHER']

    assert count_families(
        pedigree=quad_ped,
        samples=ped_samples,
    ) == {'affected': 1, 'male': 3, 'female': 1, 'quads': 1}
