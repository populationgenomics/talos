"""
script testing methods within reanalysis/validate_categories.py
"""

from dataclasses import dataclass

from reanalysis.validate_categories import gene_clean_results, count_families


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
    var_data: PicoVariant


def test_gene_clean_results():
    """
    tests the per-participant gene-filtering of results
    messy test, write and pass file paths
    """
    cat_1 = PicoVariant(['1'])
    cat_2 = PicoVariant(['2'])

    dirty_data = {
        'sam1': [PicoReport('ENSG1', cat_1), PicoReport('ENSG2', cat_1)],
        'sam2': [PicoReport('ENSG3', cat_1), PicoReport('ENSG3', cat_1)],
        'sam3': [PicoReport('ENSG4', cat_2), PicoReport('ENSG5', cat_1)],
    }
    panel_genes = {
        'genes': {
            'ENSG1': {'panels': ['1'], 'new': []},
            'ENSG2': {'panels': [], 'new': []},
            'ENSG3': {'panels': ['2'], 'new': []},
            'ENSG4': {'panels': ['3'], 'new': ['3']},
            'ENSG5': {'panels': ['3'], 'new': []},
        }
    }

    party_panels = {
        'sam1': {'panels': ['1']},
        'sam2': {'panels': ['2']},
        'sam3': {'panels': ['3']},
    }

    clean = gene_clean_results(party_panels, panel_genes, dirty_data)
    assert len(clean['sam1']) == 1
    assert clean['sam1'][0].gene == 'ENSG1'
    assert len(clean['sam2']) == 2
    assert clean['sam2'][0].gene == 'ENSG3'
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
