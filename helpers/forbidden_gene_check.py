#!/usr/bin/env python3


"""
date based forbidden gene finding
"""

import json
import logging
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

import click
import requests

from cpg_utils import to_path

PANELAPP_BASE = 'https://panelapp.agha.umccr.org/api/v1/panels'


def get_json_response(url: str) -> Any:
    """
    takes a request URL, checks for healthy response, returns the JSON
    For this purpose we only expect a dictionary return
    List use-case (activities endpoint) no longer supported

    Args:
        url (): str URL to retrieve JSON format data from

    Returns:
        the JSON response from the endpoint
    """

    response = requests.get(url, headers={'Accept': 'application/json'}, timeout=60)
    response.raise_for_status()
    return response.json()


def read_json_from_path(bucket_path: str | Path | None, default: Any = None) -> Any:
    """
    take a path to a JSON file, read into an object
    if the path doesn't exist - return the default object

    Args:
        bucket_path (str):
        default (Any):

    Returns:
        either the object from the JSON file, or None
    """

    if bucket_path is None:
        return default

    if isinstance(bucket_path, str):
        bucket_path = to_path(bucket_path)

    if bucket_path.exists():
        with bucket_path.open() as handle:
            return json.load(handle)
    return default


def get_panel_green(panel_id: int, version: str | None = None) -> set[str]:
    """
    Takes a panel number, and pulls all GRCh38 gene details from PanelApp
    For each gene, keep the MOI, symbol, ENSG (where present)

    Args:
        panel_id (): specific panel to query for
        version (): version, optional. Latest panel unless stated
    """

    # include the version if required
    panel_url = f'{PANELAPP_BASE}/{panel_id}' + (
        f'?version={version}' if version else ''
    )

    # iterate over the genes in this panel result
    return {
        gene.get('entity_name')
        for gene in get_json_response(panel_url)['genes']
        if (gene['confidence_level'] == '3' and gene['entity_type'] == 'gene')
    }


def read_panels_from_participant_file(panel_json: str) -> set[int]:
    """
    reads the per-participants panels into a set
    Args:
        panel_json (): path to a per-participant panel dump

    Returns:
        set of all the panels across all participants
    """
    participant_panels = read_json_from_path(panel_json)
    panel_set = set()
    for details in participant_panels.values():
        panel_set.update(details.get('panels', []))

    return panel_set


def find_version(panel_id: int, all_dates: list[str]) -> dict[str, str | None]:
    """
    take a panel ID and multiple date thresholds
    return the latest version of this panel at that date
    return None if no version found
    """
    panel_versions = {}

    # query for data from this endpoint
    activities: list[dict] = get_json_response(f'{PANELAPP_BASE}/{panel_id}/activities')

    for date_string in all_dates:
        date_done = False
        date_threshold = datetime.strptime(date_string, '%Y-%m-%d')

        # iterate through all activities on this panel
        for activity in activities:

            # cast the activity datetime to day-resolution
            activity_date = datetime.strptime(
                activity['created'].split('T')[0], '%Y-%m-%d'
            )

            # keep going until we land on the day, or skip past it
            if activity_date > date_threshold:
                continue

            panel_versions[date_string] = activity['panel_version']
            date_done = True
            break

        # it's possible we won't find one for some panels
        if not date_done:
            panel_versions[date_string] = None

    return panel_versions


@click.command()
@click.option('--panels', help='JSON of per-participant panels (optional)')
@click.option('--out_path', required=True, help='destination for results (directory)')
@click.argument('dates', nargs=-1)  # "YYYY-MM-DD, can be multiple whitespace-delimited"
def main(panels: str | None, out_path: str, dates: list[str]):
    """
    This script takes one or more dates, and optionally some pheno-matched data
    We first aggregate all the panels to check across all participants

    Apologies in advance to the PanelApp API - many queries

    1. Query for all panels to get current union of gene IDs across all panels

    2. For each date in turn, find the latest version of those component panels
        at that date

    3. Query for those specific panel versions, get the union, then get diff vs. new

    4. Write that result to a file

    Args:
        panels ():
        out_path ():
        dates ():
    """

    logging.info('Starting PanelApp Query Stage')
    assert dates, 'whats the point doing this with no dates?'

    if panels:
        all_panels = read_panels_from_participant_file(panels)
    else:
        all_panels = {137}

    logging.info(f'Panels: {all_panels}')

    latest_genes = set()
    for panel in all_panels:
        latest_genes.update(get_panel_green(panel))

    logging.info(f'total current genes: {len(latest_genes)}')

    # this returns {panel_id: {date: version } }
    panel_versions = {
        panel_id: find_version(panel_id, dates) for panel_id in all_panels
    }

    # check over all dates
    for date in dates:

        logging.info(f'Running the date {date}')

        date_genes = set()
        for panel, versions in panel_versions.items():
            if panel_version := versions.get(date):
                date_genes.update(get_panel_green(panel, version=panel_version))

        date_diff = sorted(latest_genes - date_genes)
        logging.info(f'date-forbidden at {date}: {len(date_diff)}')

        filename = os.path.join(out_path, f'{date}_forbidden.json')
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        with open(filename, 'w', encoding='utf-8') as handle:
            json.dump(date_diff, handle, indent=True)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )

    main()  # pylint: disable=no-value-for-parameter
