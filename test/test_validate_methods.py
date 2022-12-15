"""
script testing methods within reanalysis/validate_categories.py
"""

from dataclasses import dataclass, field

from reanalysis.validate_categories import clean_and_filter, count_families


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


def test_gene_clean_results_no_personal():
    """
    tests the per-participant gene-filtering of results
    messy test, write and pass file paths
    """

    clean = clean_and_filter(dirty_data, panel_genes, None)
    assert len(clean['sam1']) == 1
    assert clean['sam1'][0].gene == 'ENSG1'
    assert clean['sam1'][0].flags == []
    assert 'sam2' not in clean
    assert len(clean['sam3']) == 2
    assert {x.gene for x in clean['sam3']} == {'ENSG4', 'ENSG5'}


def test_gene_clean_results_personal():
    """
    tests the per-participant gene-filtering of results
    messy test, write and pass file paths
    """

    personal_panels = {
        'sam1': {'panels': [1]},
        'sam2': {'panels': []},
        'sam3': {'panels': [3, 4]},
    }

    clean = clean_and_filter(dirty_data, panel_genes, personal_panels)
    assert len(clean['sam1']) == 1
    assert clean['sam1'][0].gene == 'ENSG1'
    assert clean['sam1'][0].flags == [1]
    assert 'sam2' not in clean
    assert len(clean['sam3']) == 2
    for event in clean['sam3']:
        if event.gene == 'ENSG4':
            assert event.flags == [3]
        if event.gene == 'ENSG5':
            assert event.flags == [4]


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
        '3': 2,
    }
