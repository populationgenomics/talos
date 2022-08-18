#!/usr/bin/env python3


"""
PanelApp Forbidden gene parser
Targets the Incidentalome, and finds all genes
- green
- not tagged as cardiac

Makes them into a JSON list of forbidden content
(to be scrubbed out of any report)
"""

import json
import sys

import requests


if __name__ == '__main__':
    inci = requests.get('https://panelapp.agha.umccr.org/api/v1/panels/126')
    inci.raise_for_status()
    panel_json = inci.json()

    problem_genes = []
    for gene in panel_json['genes']:
        if gene['confidence_level'] != '3' or gene['entity_type'] != 'gene':
            continue
        if gene['tags'] != ['cardiac']:
            for build, content in gene['gene_data']['ensembl_genes'].items():
                if build.lower() == 'grch38':
                    # the ensembl version may alter over time, but will be singular
                    ensg = content[list(content.keys())[0]]['ensembl_id']
                    if not ensg:
                        print(f'FAIL! {gene}')
                        sys.exit(1)
                    problem_genes.append(ensg)

    with open('problem_genes.json', 'w', encoding='utf-8') as handle:
        json.dump(sorted(problem_genes), handle, indent=True)
