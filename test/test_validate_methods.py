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


def test_gene_clean_results_no_personal():
    """
    tests the per-participant gene-filtering of results
    messy test, write and pass file paths
    """
    cat_1 = PicoVariant(['1'])
    cat_2 = PicoVariant(['2'])
    cat_none = PicoVariant([])

    dirty_data = [
        PicoReport('ENSG1', 'sam1', cat_none),
        PicoReport('ENSG2', 'sam1', cat_1),
        PicoReport('ENSG3', 'sam2', cat_none),
        PicoReport('ENSG4', 'sam3', cat_2),
        PicoReport('ENSG5', 'sam3', cat_1),
    ]
    panel_genes = {
        'ENSG1': {'panels': [137, 1], 'new': []},
        'ENSG2': {'panels': [], 'new': []},
        'ENSG3': {'panels': [2], 'new': []},
        'ENSG4': {'panels': [3], 'new': [3]},
        'ENSG5': {'panels': [4], 'new': []},
    }

    clean = clean_and_filter(dirty_data, panel_genes, None)
    assert len(clean['sam1']) == 1
    assert {x.gene for x in clean['sam1']} == {'ENSG2'}
    assert 'sam2' not in clean
    assert len(clean['sam3']) == 2
    assert {x.gene for x in clean['sam3']} == {'ENSG4', 'ENSG5'}


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
