#!/usr/bin/env python3


"""
date based forbidden gene finding
"""


import json
import logging
import os
import sys
from datetime import datetime

import click

from reanalysis.utils import get_json_response, read_json_from_path


PANELAPP_BASE = 'https://panelapp.agha.umccr.org/api/v1/panels'


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


def find_version(panel_id: int, date_threshold: str) -> str | None:
    """
    take a panel ID and date threshold
    return the latest version of this panel at that date
    """

    date_threshold = datetime.strptime(date_threshold, '%Y-%m-%d')

    # query for data from this endpoint
    activities: list[dict] = get_json_response(f'{PANELAPP_BASE}/{panel_id}/activities')

    # iterate through all activities on this panel
    for activity in activities:

        # cast the activity datetime to day-resolution
        activity_date = datetime.strptime(activity['created'].split('T')[0], '%Y-%m-%d')

        # keep going until we land on the day, or skip past it
        if activity_date > date_threshold:
            continue

        return activity['panel_version']

    # it's possible we won't find one for some panels
    return None


@click.command()
@click.option('--panels', help='JSON of per-participant panels (optional)')
@click.option('--outpath', required=True, help='destination for results (directory)')
@click.argument('dates', nargs=-1)  # "YYYY-MM-DD, can be multiple whitespace-delimited"
def main(panels: str | None, outpath: str, dates: list[str]):
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
        outpath ():
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

    # ok, breathe...
    for date in dates:

        logging.info(f'Running the date {date}')

        panel_versions = {
            panel_id: find_version(panel_id, date) for panel_id in all_panels
        }

        date_genes = set()
        for panel, version in panel_versions.items():
            if version is None:
                continue
            date_genes.update(get_panel_green(panel, version=version))

        date_diff = sorted(latest_genes - date_genes)
        logging.info(f'date-forbidden at {date}: {len(date_diff)}')

        with open(
            os.path.join(outpath, f'{date}_forbidden.json'), 'w', encoding='utf-8'
        ) as handle:
            json.dump(date_diff, handle, indent=True)


if __name__ == '__main__':

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )

    main()  # pylint: disable=no-value-for-parameter
