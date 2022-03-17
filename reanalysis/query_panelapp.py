#!/usr/bin/env python3


"""
PanelApp Parser for Reanalysis project

Takes a panel ID
For the latest content, pulls Symbol, ENSG, and MOI

Optionally user can provide a panel version number in the past
Pull all details from the earlier version
Store all discrepancies between earlier and current

Write all output to a JSON dictionary
"""


from typing import Any, Dict, List, Optional, Union, Set
import logging
import json
import requests

from cloudpathlib import CloudPath

import click

# panelapp URL constants
PANEL_CONTENT = 'https://panelapp.agha.umccr.org/api/v1/{panel_id}'


def parse_gene_list(path_to_list: str) -> Set[str]:
    """
    parses a file (GCP or local), extracting a set of genes
    required format: clean data, one per line, gene name only
    :param path_to_list:
    :return:
    """
    if path_to_list.startswith('gs://'):
        # handle as a cloud file
        handle = CloudPath(path_to_list)

    else:

        handle = open(path_to_list, 'r', encoding='utf-8')  # pylint: disable=R1732

    return {line.rstrip() for line in handle if line.rstrip() != ''}


def get_json_response(url: str) -> Union[List[Dict[str, str]], Dict[str, Any]]:
    """
    takes a request URL, checks for healthy response, returns the JSON
    :param url:
    :return:
    """

    response = requests.get(url, headers={'Accept': 'application/json'})
    response.raise_for_status()
    return response.json()


def get_panel_green(
    panel_id: str = '137',
    version: Optional[str] = None,
) -> Dict[str, Dict[str, Union[str, bool]]]:
    """
    Takes a panel number, and pulls all GRCh38 gene details from PanelApp
    For each gene, keep the MOI, symbol, ENSG (where present)

    :param panel_id: defaults to the PanelAppAU Mendeliome
    :param version:
    """

    # prepare the query URL
    panel_app_genes_url = PANEL_CONTENT.format(panel_id=panel_id)
    if version is not None:
        panel_app_genes_url += f'?version={version}'

    panel_response = requests.get(panel_app_genes_url)
    panel_response.raise_for_status()
    panel_json = panel_response.json()

    panel_version = panel_json.get('version')

    # pop the version in logging if not manually defined
    if version is None:
        logging.info('Current panel version: %s', panel_version)

    gene_dict = {'panel_metadata': {'current_version': panel_version}}

    for gene in panel_json['genes']:

        # only retain green genes
        if gene['confidence_level'] != '3' or gene['entity_type'] != 'gene':
            continue

        ensg = None
        symbol = gene.get('entity_name')

        # take the PanelApp MOI, don't simplify
        moi = gene.get('mode_of_inheritance', None)

        # for some reason the build is capitalised oddly in panelapp
        # at least one entry doesn't have an ENSG annotation
        for build, content in gene['gene_data']['ensembl_genes'].items():
            if build.lower() == 'grch38':
                # the ensembl version may alter over time, but will be singular
                ensg = content[list(content.keys())[0]]['ensembl_id']

        # this appears to be missing in latest panel version
        if symbol == 'RNU12' and ensg is None:
            ensg = 'ENSG00000276027'

        # save the entity into the final dictionary
        # include fields to recognise altered gene data
        gene_dict[ensg] = {
            'symbol': symbol,
            'moi': moi,
            'new': False,
            'changed': False,
            'old_moi': None,
        }

    return gene_dict


def get_panel_changes(
    previous_version: str,
    panel_id: str,
    latest_content: Dict[str, Dict[str, Union[str, bool]]],
):
    """
    take the latest panel content, and compare with a previous version
    update content in original dict where appropriate
    https://panelapp.agha.umccr.org/api/v1/panels/137/?version=0.10952

    :param previous_version:
    :param panel_id:
    :param latest_content:
    :return: None, updates object in place
    """

    # get the full content for the specified panel version
    previous_content = get_panel_green(panel_id=panel_id, version=previous_version)

    # iterate over the latest content,skipping over the metadata keys
    for gene_ensg in [
        ensg for ensg in latest_content.keys() if ensg != 'panel_metadata'
    ]:

        value = latest_content[gene_ensg]

        # if the gene wasn't present before, take it in full
        if gene_ensg not in previous_content:
            latest_content[gene_ensg]['new'] = True

        # otherwise check if the MOI has changed
        else:
            prev_moi = previous_content.get(gene_ensg).get('moi')
            latest_moi = value.get('moi')

            # if so, store the old and new MOI
            if prev_moi != latest_moi:
                latest_content[gene_ensg]['changed'] = True
                latest_content[gene_ensg]['old_moi'] = prev_moi


def gene_list_differences(
    latest_content: Dict[str, Dict[str, Union[str, bool]]], previous_genes: Set[str]
):
    """
    takes a gene list representing prior data,
    identifies genes as 'new' where absent in that reference data

    :param latest_content:
    :param previous_genes:
    :return: None, updates object in place
    """
    for gene_ensg in [
        ensg for ensg in latest_content.keys() if ensg != 'panel_metadata'
    ]:
        if gene_ensg not in previous_genes:
            latest_content[gene_ensg]['new'] = True


@click.command()
@click.option('--id', 'panel_id', default='137', help='ID to use in panelapp')
@click.option('--out', 'out_path', help='path to write resulting JSON to')
@click.option(
    '--previous',
    help='identify panel differences between current version and this previous version',
)
@click.option('--gene_list', help='pointer to a file, containing a prior gene list')
def main(
    panel_id: str,
    out_path: str,
    previous: Optional[str] = None,
    gene_list: Optional[str] = None,
):
    """
    takes a panel ID
    finds all latest panel data from the API
    optionally also takes a prior version argument, records all panel differences
        - new genes
        - altered MOI
    :param panel_id:
    :param out_path: path to write a JSON object out to
    :param previous: prior panel version to compare to
    :param gene_list: alternative to prior data, give a strict gene list file
    :return:
    """

    if gene_list is not None and previous is not None:
        raise ValueError('Only one of [Date/GeneList] can be specified per run')

    # get latest panel data
    panel_dict = get_panel_green(panel_id=panel_id)

    if gene_list is not None:
        gene_list_contents = parse_gene_list(gene_list)
        gene_list_differences(panel_dict, gene_list_contents)

    # migrate more of this into a method to test
    if previous is not None:
        # only continue if the versions are different
        if previous != panel_dict['panel_metadata'].get('current_version'):
            logging.info('Previous panel version: %s', previous)
            get_panel_changes(
                previous_version=previous,
                panel_id=panel_id,
                latest_content=panel_dict,
            )
            panel_dict['panel_metadata']['previous_version'] = previous

    logging.info('Writing output JSON file to %s', out_path)
    with open(out_path, 'w', encoding='utf-8') as handle:
        json.dump(panel_dict, handle, indent=True, default=str)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=E1120
