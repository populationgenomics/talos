"""
script testing methods within reanalysis/validate_categories.py
"""

import json
import os
from dataclasses import dataclass

from reanalysis.validate_categories import gene_clean_results, update_result_meta


@dataclass
class PicoVariant:
    """
    smallest Reportable variant object for this test
    """

    gene: str


def test_gene_clean_results(tmpdir):
    """
    tests the per-participant gene-filtering of results
    messy test, write and pass file paths
    """
    dirty_data = {
        'metadata': {'foo': 'bar'},
        'sam1': [PicoVariant('ENSG1'), PicoVariant('ENSG2')],
        'sam2': [PicoVariant('ENSG3'), PicoVariant('ENSG3')],
        'sam3': [PicoVariant('ENSG4'), PicoVariant('ENSG5')],
    }
    panel_genes = os.path.join(tmpdir, 'panel_genes.json')
    with open(panel_genes, 'w', encoding='utf-8') as handle:
        json.dump({'default': ['ENSG1'], '2': ['ENSG3'], '3': ['ENSG5']}, handle)

    party_panels = os.path.join(tmpdir, 'party_panels.json')
    with open(party_panels, 'w', encoding='utf-8') as handle:
        json.dump({'sam2': {'panels': ['2']}, 'sam3': {'panels': ['3']}}, handle)

    clean = gene_clean_results(party_panels, panel_genes, dirty_data)
    assert len(clean['sam1']) == 1
    assert clean['sam1'][0].gene == 'ENSG1'
    assert len(clean['sam2']) == 2
    assert clean['sam2'][0].gene == 'ENSG3'
    assert len(clean['sam3']) == 1
    assert clean['sam3'][0].gene == 'ENSG5'


def test_update_results_meta(peddy_ped):
    """
    testing the dict update
    """

    results = {}
    panelapp = {
        'metadata': [
            {'name': 'biff', 'version': 'pow', 'id': 'wallop'},
            {'name': 'extra_panel', 'version': 'extra_version', 'id': 2},
        ]
    }

    ped_samples = ['male', 'female', 'mother_1', 'father_1', 'mother_2', 'father_2']

    big_results = update_result_meta(
        results=results,
        pedigree=peddy_ped,
        panelapp=panelapp,
        samples=ped_samples,
    )
    assert big_results == {
        'metadata': {
            'run_datetime': 'bar',
            'input_file': 'foo',
            'cohort': 'cohort',
            'family_breakdown': {
                'affected': 2,
                'male': 3,
                'female': 3,
                'trios': 2,
                '3': 2,
            },
            'panels': [
                {'name': 'biff', 'version': 'pow', 'id': 'wallop'},
                {'name': 'extra_panel', 'version': 'extra_version', 'id': 2},
            ],
        }
    }
