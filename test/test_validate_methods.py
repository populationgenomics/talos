"""
script testing methods within reanalysis/validate_categories.py
"""


from dataclasses import dataclass, field

from reanalysis.validate_categories import (
    clean_and_filter,
    count_families,
    prepare_results_shell,
)


@dataclass
class PicoVariant:
    """
    some categories
    """

    categories: list


@dataclass
class PicoReport:
    """
    smallest Reportable variant object for this test
    """

    gene: str
    sample: str
    var_data: PicoVariant
    reasons: set[str] = field(default_factory=set)
    flags: list[str] = field(default_factory=list)
    panels: dict[str] = field(default_factory=dict)

    def __lt__(self, other):
        return True

    def __eq__(self, other):
        return False


cat_1 = PicoVariant(['1'])
cat_2 = PicoVariant(['2'])
cat_none = PicoVariant([])

dirty_data = [
    PicoReport('ENSG1', 'sam1', cat_1),
    PicoReport('ENSG2', 'sam1', cat_none),
    PicoReport('ENSG3', 'sam2', cat_none),
    PicoReport('ENSG4', 'sam3', cat_2),
    PicoReport('ENSG5', 'sam3', cat_1),
]
panel_genes = {
    'metadata': [
        {'id': 137, 'name': 137},
        {'id': 1, 'name': 1},
        {'id': 2, 'name': 2},
        {'id': 3, 'name': 3},
        {'id': 4, 'name': 4},
    ],
    'genes': {
        'ENSG1': {'panels': [137, 1], 'new': []},
        'ENSG2': {'panels': [], 'new': []},
        'ENSG3': {'panels': [2], 'new': []},
        'ENSG4': {'panels': [3], 'new': [3]},
        'ENSG5': {'panels': [4], 'new': []},
    },
}


def test_results_shell(peddy_ped):
    """

    Returns:

    """
    samples = ['male', 'female', 'irrelevant']
    sample_panels = {
        'male': {'panels': [1, 3], 'hpo_terms': ['Boneitis']},
        'female': {'panels': [1, 2], 'hpo_terms': ['HPfemale']},
        'other': [4],
    }
    panelapp = {
        'metadata': [
            {'id': 1, 'name': 'lorem'},
            {'id': 2, 'name': 'ipsum'},
            {'id': 3, 'name': 'etc'},
        ]
    }
    shell = prepare_results_shell(samples, peddy_ped, sample_panels, panelapp)

    # top level only has the two affected participants
    expected = {
        'male': {
            'variants': [],
            'metadata': {
                'ext_id': 'male',
                'family_id': 'family_1',
                'members': {
                    'male': {'sex': 'male', 'affected': True, 'ext_id': 'male'},
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
            'variants': [],
            'metadata': {
                'ext_id': 'female',
                'family_id': 'family_2',
                'members': {
                    'female': {'sex': 'female', 'affected': True, 'ext_id': 'female'},
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
            },
        },
    }

    assert shell == expected


def test_gene_clean_results_no_personal():
    """
    tests the per-participant gene-filtering of results
    messy test, write and pass file paths
    """
    results_holder = {
        'sam1': {'variants': []},
        'sam2': {'variants': []},
        'sam3': {'variants': []},
    }

    clean = clean_and_filter(
        results_holder=results_holder, result_list=dirty_data, panelapp_data=panel_genes
    )
    assert len(clean['sam1']['variants']) == 1
    assert clean['sam1']['variants'][0].gene == 'ENSG1'
    assert clean['sam1']['variants'][0].flags == []
    assert len(clean['sam2']['variants']) == 0
    assert len(clean['sam3']['variants']) == 2
    assert {x.gene for x in clean['sam3']['variants']} == {'ENSG4', 'ENSG5'}


def test_gene_clean_results_personal():
    """
    tests the per-participant gene-filtering of results
    messy test, write and pass file paths
    """
    results_holder = {
        'sam1': {'variants': []},
        'sam2': {'variants': []},
        'sam3': {'variants': []},
    }
    personal_panels = {
        'sam1': {'panels': [1], 'hpo_terms': ['HP1']},
        'sam2': {'panels': [], 'hpo_terms': ['HP2']},
        'sam3': {'panels': [3, 4], 'hpo_terms': ['HP3']},
    }

    clean = clean_and_filter(results_holder, dirty_data, panel_genes, personal_panels)
    assert len(clean['sam1']['variants']) == 1
    assert clean['sam1']['variants'][0].gene == 'ENSG1'
    assert len(clean['sam1']['variants'][0].flags) == 0
    assert clean['sam1']['variants'][0].panels['matched'] == [1]
    assert len(clean['sam2']['variants']) == 0
    assert len(clean['sam3']['variants']) == 2
    for event in clean['sam3']['variants']:
        if event.gene == 'ENSG4':
            assert event.panels['matched'] == [3]
        if event.gene == 'ENSG5':
            assert event.panels['matched'] == [4]


def test_update_results_meta(peddy_ped):
    """
    testing the dict update
    """

    ped_samples = ['male', 'female', 'mother_1', 'father_1', 'mother_2', 'father_2']

    assert count_families(pedigree=peddy_ped, samples=ped_samples,) == {
        'affected': 2,
        'male': 3,
        'female': 3,
        'trios': 2,
    }


def test_update_results_missing_father(peddy_ped):
    """
    testing the dict update
    """

    ped_samples = ['male', 'female', 'mother_1', 'mother_2', 'father_2']

    assert count_families(
        pedigree=peddy_ped,
        samples=ped_samples,
    ) == {'affected': 2, 'male': 2, 'female': 3, 'trios': 1, '2': 1}


def test_update_results_quad(quad_ped):
    """
    testing the dict update
    """

    ped_samples = ['PROBAND', 'SIBLING', 'FATHER', 'MOTHER']

    assert count_families(
        pedigree=quad_ped,
        samples=ped_samples,
    ) == {'affected': 1, 'male': 3, 'female': 1, 'quads': 1}
