#!/usr/bin/env python3


"""
PanelApp Parser for Reanalysis project

 Takes a panel ID
Pulls latest 'green' content; Symbol, ENSG, and MOI

Optionally user can provide a panel version number in the past
Pull all details from the earlier version
Annotate all discrepancies between earlier and current

Optionally user can provide path to a JSON gene list
Annotate all genes in current panel and not the gene list

Write all output to a JSON dictionary
"""


from typing import Any, Union
import logging
import json
import sys

from argparse import ArgumentParser

from cpg_utils import to_path

from reanalysis.utils import get_json_response, read_json_from_path


MENDELIOME = '137'
PANELAPP_BASE = 'https://panelapp.agha.umccr.org/api/v1/panels/'
PanelData = dict[str, dict | list[dict]]


def parse_gene_list(path_to_list: str) -> set[str]:
    """
    parses a json file (GCP or local), extracting a set of genes
    required format: a json list of strings
    :param path_to_list:
    """
    logging.info(f'Loading gene list from {path_to_list}')
    with open(to_path(path_to_list), encoding='utf-8') as handle:
        return set(json.load(handle))


def get_panel_green(
    panel_id: str = MENDELIOME,
) -> dict[str, dict[str, Union[str, bool]]]:
    """
    Takes a panel number, and pulls all GRCh38 gene details from PanelApp
    For each gene, keep the MOI, symbol, ENSG (where present)

    :param panel_id: defaults to the PanelAppAU Mendeliome
    """

    # prepare the query URL
    panel_app_genes_url = f'{PANELAPP_BASE}{panel_id}'
    panel_json = get_json_response(panel_app_genes_url)

    panel_version = panel_json.get('version')
    panel_name = panel_json.get('name')

    # pop the version in logging if not manually defined
    logging.info(
        f'Current {panel_name} panel version: {panel_version}',
    )

    gene_dict = {
        'metadata': [
            {
                'name': panel_name,
                'version': panel_version,
                'id': panel_id,
            }
        ]
    }

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

        if ensg is None:
            logging.info(f'Gene "{symbol} lacks an ENSG ID, so it is being excluded')
            continue

        # save the entity into the final dictionary
        # include fields to recognise altered gene data
        gene_dict[ensg] = {
            'symbol': symbol,
            'moi': moi,
            'new': False,
            'flags': [],
        }

    return gene_dict


def gene_list_differences(latest_content: PanelData, previous_genes: set[str]):
    """
    takes a gene list representing prior data,
    identifies genes as 'new' where absent in that reference data

    :param latest_content:
    :param previous_genes:
    :return: None, updates object in place
    """
    new_genes = 0
    for content in [
        content for ensg, content in latest_content.items() if ensg != 'metadata'
    ]:
        if content['symbol'] not in previous_genes:
            content['new'] = True
            new_genes += 1
    logging.info(f'{new_genes} genes were not present in the gene list')


def write_output_json(output_path: str, object_to_write: Any):
    """
    writes object to a json file
    to_path provides platform abstraction

    :param output_path:
    :param object_to_write:
    """

    logging.info(f'Writing output JSON file to {output_path}')
    out_route = to_path(output_path)

    if out_route.exists():
        logging.info(f'Output path "{output_path}" exists, will be overwritten')

    with out_route.open('w') as fh:
        json.dump(object_to_write, fh, indent=4, default=str)


def combine_mendeliome_with_other_panels(panel_dict: PanelData, additional: PanelData):
    """
    takes the main panel data and an additional panel dict

    :param panel_dict:
    :param additional:
    """

    additional_name = additional['metadata'][0]['name']
    panel_dict['metadata'].append(
        {
            'name': additional['metadata'][0]['name'],
            'version': additional['metadata'][0]['version'],
            'id': additional['metadata'][0]['id'],
        }
    )
    panel_keys = panel_dict.keys()
    for ensg in additional.keys():
        if ensg == 'metadata':
            continue

        if ensg in panel_keys:
            # update MOI if None
            if panel_dict[ensg]['moi'] is None:
                panel_dict[ensg]['moi'] = additional[ensg].get('moi', None)
            panel_dict[ensg]['flags'].append(additional_name)

        else:
            panel_dict[ensg] = {
                'symbol': additional[ensg].get('symbol'),
                'moi': additional[ensg].get('moi', None),
                'new': False,
                'flags': [additional_name],
            }


def get_list_from_participants(participant_panels) -> set[str]:
    """
    takes the per-participant panels and extracts the panel list

    Parameters
    ----------
    participant_panels : the dictionary of per-participant fam_id, hpo & panels

    Returns
    -------
    list of all unique panels
    """

    panel_set = set()
    for details in participant_panels.values():
        panel_set.update(details.get('panels', []))

    return panel_set


def grab_genes_only(panel_data: PanelData) -> list[str]:
    """
    takes panelapp data for a panel and grabs the gene IDs

    Parameters
    ----------
    panel_data : a full set of data from a panel

    Returns
    -------
    a list of only the ENSG IDs
    """

    return [key for key in panel_data.keys() if key != 'metadata']


def main(
    panel_list: list[str] | set[str],
    panel_file: str,
    out_path: str,
    gene_list: str | None,
):
    """
    Base assumption here is that we are always using the Mendeliome
    Optionally, additional panel IDs can be specified to expand the gene list

    Finds all latest panel data from the API
    optionally take a prior version argument, records all panel differences
        - new genes
        - altered MOI
    optionally take a reference to a JSON gene list, records all genes:
        - green in current panelapp
        - absent in provided gene list
    :param panel_list: iterable of panelapp IDs
    :param panel_file: json file of panel IDs per participant
    :param out_path: path to write a JSON object out to
    :param gene_list: alternative to prior data, give a strict gene list file
    :return:
    """

    logging.info('Starting PanelApp Query Stage')

    # get latest Mendeliome data
    panel_dict = get_panel_green()

    # create this to stop the linter complaining
    panel_master = {}

    # we can merge functionality here
    if panel_file:
        # load the json dictionary
        participant_panels = read_json_from_path(panel_file)
        panel_list = get_list_from_participants(participant_panels)
        panel_master = {'default': grab_genes_only(panel_dict)}

    if panel_list:
        for additional_panel_id in panel_list:
            ad_panel = get_panel_green(panel_id=additional_panel_id)
            combine_mendeliome_with_other_panels(
                panel_dict=panel_dict, additional=ad_panel
            )

            # if we are finding per-panel gene lists, store the genes for
            # this panel in the X-cohort master dict
            if panel_file:
                panel_master[str(additional_panel_id)] = grab_genes_only(ad_panel)

    if gene_list is not None:
        logging.info(f'A Gene_List was selected: {gene_list}')
        gene_list_contents = parse_gene_list(gene_list)
        logging.info(f'Length of gene list: {len(gene_list_contents)}')
        gene_list_differences(panel_dict, gene_list_contents)

    write_output_json(f'{out_path}.json', panel_dict)

    # if we populated panel master content - write that too
    # trying not to make the API too messy... but this is an optional output file
    if panel_master:
        write_output_json(f'{out_path}_per_panel.json', panel_master)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )

    parser = ArgumentParser()
    panel_input = parser.add_mutually_exclusive_group()
    panel_input.add_argument_group()
    panel_input.add_argument(
        '-p',
        nargs='+',
        default=[],
        help='Panelapp IDs of any additional panels to query for',
    )
    panel_input.add_argument(
        '--panel_file', type=str, help='json file containing per-participant panels'
    )

    parser.add_argument(
        '--out_path', type=str, required=True, help='Path to write output JSON to'
    )
    parser.add_argument(
        '--gene_list',
        type=str,
        required=False,
        help='If a gene list is being used as a comparison',
    )
    args = parser.parse_args()
    main(
        panel_list=args.p,
        panel_file=args.panel_file,
        out_path=args.out_path,
        gene_list=args.gene_list,
    )
