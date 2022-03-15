#!/usr/bin/env python3


"""
PanelApp Parser for Reanalysis

Takes a panel ID
For the latest content, pulls Symbol, ENSG, and MOI
    (MOI is simplified from PanelApp enum)

Optionally user can provide a date in the past
Identify the highest panel version prior to that date
Pull all details from the earlier version
Store all discrepancies between earlier and current

Write all output to a JSON dictionary
"""


import logging
from datetime import datetime
from typing import Any, Dict, List, Optional, Union
import json

import requests

import click


# panelapp URL constants
PANELAPP_ROOT = 'https://panelapp.agha.umccr.org/api/v1'
PANEL_ROOT = f'{PANELAPP_ROOT}/panels'
PANEL_CONTENT = f'{PANEL_ROOT}/{{panel_id}}'
ACTIVITIES = f'{PANEL_CONTENT}/activities'


def get_json_response(url: str) -> Union[List[Dict[str, str]], Dict[str, Any]]:
    """
    takes a request URL, checks for healthy response, returns the JSON
    :param url:
    :return:
    """

    response = requests.get(url, headers={'Accept': 'application/json'})
    response.raise_for_status()
    return response.json()


def get_previous_version(panel_id: str, since: datetime) -> str:
    """
    work through the list of panel updates in reverse (earliest -> latest)
    return the first panel version after the threshold date
    If we don't find any, return the latest version

    Note: Django dates include Hours and Minutes, so any Django date from
    the same day as a 'YYYY-MM-DD' date will be 'greater than'
    for the purposes of a value comparison

    revisit the exact implementation
    consider replacement with a panel version (instead of date)

    :param panel_id: panel ID to use
    :param since: date of the
    :return:
    """

    activity_list = get_json_response(url=ACTIVITIES.format(panel_id=panel_id))
    entry_version = None
    for entry in reversed(activity_list):
        # take the panel version
        entry_version = entry.get('panel_version')
        # uses the django datestamp format
        entry_date = datetime.strptime(entry.get('created'), '%Y-%m-%dT%H:%M:%S.%fz')

        if entry_date >= since:
            return entry_version
    return entry_version


def get_panel_green(
    reference_genome: Optional[str] = 'GRch38',
    panel_id: str = '137',
    version: Optional[str] = None,
) -> Dict[str, Dict[str, Union[str, bool]]]:
    """
    Takes a panel number, and pulls all details from PanelApp
    :param reference_genome: GRch37 or GRch38
    :param panel_id: defaults to the PanelAppAU Mendeliome
    :param version:
    :return:
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
            if build.lower() == reference_genome.lower():
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
    :return:
    """

    # get the full content for the specified panel version
    previous_content = get_panel_green(panel_id=panel_id, version=previous_version)
    # iterate over the latest content
    # skip over the metadata keys
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


@click.command()
@click.option('--id', 'panel_id', default='137', help='ID to use in panelapp')
@click.option('--out', 'out_path', help='path to write resulting JSON to')
@click.option(
    '--date',
    help='identify panel differences between this date and now (YYYY-MM-DD)',
)
def main(panel_id: str, out_path: str, date: Optional[str] = None):
    """
    takes a panel ID and a date
    finds all latest panel data from the API
    uses activities endpoint for highest panel version prior to _date_
    retrieves panel data at that version
    compares, and records all panel differences
        - new genes
        - altered MOI
    :param panel_id:
    :param out_path: path to write a JSON object out to
    :param date: string to parse as a Datetime
    :return:
    """

    # get latest panel data
    panel_dict = get_panel_green(panel_id=panel_id)

    # migrate more of this into a method to test
    if date is not None:
        since_datetime = datetime.strptime(date, '%Y-%m-%d')
        if since_datetime > datetime.today():
            raise ValueError(f'The specified date {date} cannot be in the future')

        early_version = get_previous_version(panel_id=panel_id, since=since_datetime)

        # only continue if the versions are different
        if early_version != panel_dict['panel_metadata'].get('current_version'):
            logging.info('Previous panel version: %s', early_version)
            logging.info('Previous version date: %s', date)
            get_panel_changes(
                previous_version=early_version,
                panel_id=panel_id,
                latest_content=panel_dict,
            )
            panel_dict['panel_metadata']['previous_version'] = early_version
            panel_dict['panel_metadata']['previous_date'] = date

    logging.info('Writing output JSON file to %s', out_path)
    with open(out_path, 'w', encoding='utf-8') as handle:
        json.dump(panel_dict, handle, indent=True, default=str)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=E1120
