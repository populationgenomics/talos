#!/usr/bin/env python3

"""
Start over...
"""


import logging
import sys

import click

from cpg_utils.config import get_config

from reanalysis.utils import (
    read_json_from_path,
    get_json_response,
    get_simple_moi,
    write_output_json,
    ORDERED_MOIS,
)


PanelData = dict[str, dict | list[dict]]
OLD_DATA = {'genes': {}}


def is_this_gene_new(gene_id: str, panel_id: int) -> bool:
    """
    check against the historic data to see if this gene is new
    this is an interesting one - should we permit an 'all' panelID
    Args:
        gene_id (): ENSG***
        panel_id (): integer panelID

    Returns:
        True if the gene was not seen in the previous data
        otherwise, False
    """
    if gene_id not in OLD_DATA['genes']:
        return True

    # it's new if the panel ID or 'all' wasn't already associated
    return not any(
        {
            pan in OLD_DATA['genes'][gene_id].get('panels', [])
            for pan in [panel_id, 'all']
        }
    )


def request_panel_data(url: str) -> tuple[str, str, list]:
    """
    just takes care of the panelapp query
    Args:
        url ():

    Returns:
        components of the panelapp response
    """

    panel_json = get_json_response(url)

    # steal attributes
    panel_name = panel_json.get('name')
    panel_version = panel_json.get('version')
    panel_genes = panel_json.get('genes')

    # log name and version
    logging.info(f'Current {panel_name} version: {panel_version}')

    return panel_name, panel_version, panel_genes


def get_panel_green(gene_dict: PanelData = None, panel_id: int | None = None):
    """
    Takes a panel number, and pulls all GRCh38 gene details from PanelApp
    For each gene, keep the MOI, symbol, ENSG (where present)

    Args:
        gene_dict (): dictionary to continue populating
        panel_id (): specific panel or 'base' (e.g. 137)
    """

    if panel_id is None:
        panel_id = get_config()['workflow'].get('default_panel', 137)

    # prepare the query URL
    panelapp_base = get_config()['workflow'].get(
        'panelapp', 'https://panelapp.agha.umccr.org/api/v1/panels/'
    )

    panel_name, panel_version, panel_genes = request_panel_data(
        f'{panelapp_base}{panel_id}'
    )

    # add metadata for this panel & version
    gene_dict['metadata'].append(
        {
            'name': panel_name,
            'version': panel_version,
            'id': panel_id,
        }
    )

    # iterate over the genes in this panel result
    for gene in panel_genes:

        # only retain green genes
        if gene['confidence_level'] != '3' or gene['entity_type'] != 'gene':
            continue

        ensg = None
        symbol = gene.get('entity_name')

        # for some reason the build is capitalised oddly in panelapp
        # at least one entry doesn't have an ENSG annotation
        for build, content in gene['gene_data']['ensembl_genes'].items():
            if build.lower() == 'grch38':
                # the ensembl version may alter over time, but will be singular
                ensg = content[list(content.keys())[0]]['ensembl_id']

        if ensg is None:
            logging.info(f'Gene "{symbol} lacks an ENSG ID, so it is being excluded')
            continue

        # check if this is a new gene in this analysis
        new_gene = is_this_gene_new(ensg, panel_id)
        moi = get_simple_moi(gene.get('mode_of_inheritance'))
        print(ensg, new_gene, moi)

        # either update or add a new entry
        if ensg in gene_dict['genes'].keys():
            this_gene = gene_dict['genes'][ensg]

            # if we have one biallelic, and one monoallelic, merge
            if all(
                _moi in ('Biallelic', 'Monoallelic') for _moi in [moi, this_gene['moi']]
            ):
                moi = 'Mono_And_Biallelic'

            # pylint: disable=unnecessary-lambda
            # take the more lenient of the gene MOI options
            this_gene['moi'] = sorted(
                [moi, this_gene['moi']],
                key=lambda x: ORDERED_MOIS.index(x),
            )[0]

            # add this panel to the list
            this_gene['panels'].append(panel_id)

            # if this is/was new - it's new
            this_gene['new'] = this_gene['new'] or new_gene

        else:
            # save the entity into the final dictionary
            gene_dict['genes'][ensg] = {
                'symbol': symbol,
                'moi': moi,
                'new': new_gene,
                'panels': [panel_id],
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


@click.command()
@click.option('--panels', help='JSON of per-participant panels')
@click.option('--out_path', required=True, help='destination for results')
@click.option('--previous', help="previous data for finding 'new' genes")
def main(panels: str | None, out_path: str, previous: str | None):
    """
    if present, reads in any prior reference data
    if present, reads additional panels to use
    queries panelapp for each panel in turn, aggregating results

    Args:
        panels ():
        out_path ():
        previous ():
    """

    logging.info('Starting PanelApp Query Stage')

    # update OLD_DATA global object - needs to be a json in this format
    # absolutely no need for this to be global?
    # {'genes': {'ENSG***': {'panels': [1, 2, 3]}}}
    # pylint: disable=global-statement
    global OLD_DATA
    if prev := previous:
        OLD_DATA = read_json_from_path(prev)
    elif prev := get_config()['workflow'].get('pre_panelapp'):
        OLD_DATA = read_json_from_path(prev)
    else:
        logging.info('No prior data found, all genes are new...')

    # set up the gene dict
    gene_dict: PanelData = {'metadata': [], 'genes': {}}

    # first add the base content
    get_panel_green(gene_dict)

    # if participant panels were provided, add each of those to the gene data
    if panels is not None:
        panel_list = read_panels_from_participant_file(panels)
        logging.info(f'All additional panels: {", ".join(map(str, panel_list))}')
        for panel in panel_list:
            get_panel_green(gene_dict=gene_dict, panel_id=panel)

    # write the output to long term storage
    write_output_json(output_path=out_path, object_to_write=gene_dict)


if __name__ == '__main__':

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )

    main()  # pylint: disable=no-value-for-parameter
