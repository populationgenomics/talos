#!/usr/bin/env python3


"""
Complete revision
"""

# mypy: ignore-errors

import logging
import sys
from datetime import datetime
from dateutil.relativedelta import relativedelta

import click
from cpg_utils.config import get_config

from reanalysis.utils import (
    find_latest_file,
    get_cohort_config,
    get_json_response,
    get_simple_moi,
    read_json_from_path,
    save_new_historic,
    write_output_json,
    IRRELEVANT_MOI,
    ORDERED_MOIS,
)

PanelData = dict[str, dict | list[dict]]
PANELAPP_HARD_CODED_DEFAULT = 'https://panelapp.agha.umccr.org/api/v1/panels'
PANELAPP_BASE = get_config()['panels'].get('panelapp', PANELAPP_HARD_CODED_DEFAULT)
DEFAULT_PANEL = get_config()['panels'].get('default_panel', 137)


def request_panel_data(url: str) -> tuple[str, str, list]:
    """
    takes care of the panelapp query
    Args:
        url ():

    Returns:
        components of the panelapp response
    """

    panel_json = get_json_response(url)
    panel_name = panel_json.get('name')
    panel_version = panel_json.get('version')
    panel_genes = panel_json.get('genes')

    # log name and version
    logging.info(f'{panel_name} version: {panel_version}')

    return panel_name, panel_version, panel_genes


def get_panel_green(
    gene_dict: PanelData,
    old_data: dict[str, list],
    panel_id: int = DEFAULT_PANEL,
    version: str | None = None,
    blacklist: list[str] | None = None,
    forbidden_genes: set[str] | None = None,
):
    """
    Takes a panel number, and pulls all GRCh38 gene details from PanelApp
    For each gene, keep the MOI, symbol, ENSG (where present)

    Args:
        gene_dict (): dictionary to continue populating
        old_data (dict[str, list]): dict of lists - panels per gene
        panel_id (): specific panel or 'base' (e.g. 137)
        version (): version, optional. Latest panel unless stated
        blacklist (): list of symbols/ENSG IDs to remove from this panel
        forbidden_genes (set[str]): genes to remove for this cohort
    """

    if blacklist is None:
        blacklist = []
    if forbidden_genes is None:
        forbidden_genes = set()

    # include the version if required
    panel_url = f'{PANELAPP_BASE}/{panel_id}' + (
        f'?version={version}' if version else ''
    )

    panel_name, panel_version, panel_genes = request_panel_data(panel_url)

    # add metadata for this panel & version
    gene_dict['metadata'].append(
        {'name': panel_name, 'version': panel_version, 'id': panel_id}
    )

    # iterate over the genes in this panel result
    for gene in panel_genes:

        symbol = gene.get('entity_name')

        # only retain green genes
        if (
            gene['confidence_level'] != '3'
            or gene['entity_type'] != 'gene'
            or symbol in forbidden_genes
        ):
            continue

        ensg = None
        chrom = None

        # for some reason the build is capitalised oddly in panelapp
        # at least one entry doesn't have an ENSG annotation
        for build, content in gene['gene_data']['ensembl_genes'].items():
            if build.lower() == 'grch38':
                # the ensembl version may alter over time, but will be singular
                ensembl_data = content[list(content.keys())[0]]
                ensg = ensembl_data['ensembl_id']
                chrom = ensembl_data['location'].split(':')[0]

        if (
            ensg is None
            or ensg in blacklist
            or symbol in blacklist
            or ensg in forbidden_genes
        ):
            logging.info(f'Gene {symbol}/{ensg} removed from {panel_name}')
            continue

        # check if this is a new gene in this analysis
        new_gene = panel_id not in old_data.get(ensg, [])
        if new_gene:
            old_data.setdefault(ensg, []).append(panel_id)

        exact_moi = gene.get('mode_of_inheritance', 'unknown').lower()

        # either update or add a new entry
        if ensg in gene_dict['genes'].keys():
            this_gene = gene_dict['genes'][ensg]

            # now we find it on this panel
            this_gene['panels'].append(panel_id)

            # add this moi to the set
            if exact_moi not in IRRELEVANT_MOI:
                this_gene['moi'].add(exact_moi)

            # if this is/was new - it's new
            if new_gene:
                this_gene['new'].append(panel_id)

        else:

            # save the entity into the final dictionary
            gene_dict['genes'][ensg] = {
                'symbol': symbol,
                'moi': {exact_moi} if exact_moi not in IRRELEVANT_MOI else set(),
                'new': [panel_id] if new_gene else [],
                'panels': [panel_id],
                'chrom': chrom,
            }


def get_best_moi(gene_dict: dict):
    """
    From the collected set of all MOIs, take the most lenient
    If set was empty (i.e. no specific MOI found) find the default

    Default is found autosome/sex chrom aware

    Args:
        gene_dict (): the 'genes' index of the collected dict
    """

    for content in gene_dict.values():

        # accept the simplest MOI if no exact moi found
        if not content['moi']:
            content['moi'] = get_simple_moi(None, chrom=content['chrom'])
            continue

        # otherwise accept the most lenient valid MOI
        moi_set = {
            get_simple_moi(moi, chrom=content['chrom']) for moi in content['moi']
        }

        # force a combined MOI here
        if 'Biallelic' in moi_set and 'Monoallelic' in moi_set:
            content['moi'] = 'Mono_And_Biallelic'

        else:
            # take the more lenient of the gene MOI options
            content['moi'] = sorted(moi_set, key=lambda x: ORDERED_MOIS.index(x))[0]


def read_panels_from_participant_file(panel_json: str) -> set[int]:
    """
    reads the per-participants panels into a set
    Args:
        panel_json (): path to a per-participant panel dump

    Returns:
        set of all the panels across all participants
    """
    participant_panels = read_json_from_path(panel_json, {})
    panel_set = set()
    for details in participant_panels.values():
        panel_set.update(details.get('panels', []))

    return panel_set


def find_core_panel_version() -> str | None:
    """
    take the default panel ID from config
    iterate through its associated activities
    return the panel version which was closest to 12 months ago

    or return None, i.e. if the panel is not >= 12 months old

    Returns:
        a version string from X months prior - see config
    """

    date_threshold = datetime.today() - relativedelta(
        months=get_config()['panels']['panel_month_delta']
    )

    # query for data from this endpoint
    activities: list[dict] = get_json_response(
        f'{PANELAPP_BASE}/{DEFAULT_PANEL}/activities'
    )

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


def get_new_genes(
    current_genes: set[str], old_version: str, forbidden: set[str] | None = None
) -> set[str]:
    """
    query for two versions of the same panel
    find all genes new on that panel between versions
    Args:
        current_genes (): current Mendeliome genes
        old_version ():
        forbidden ():

    Returns:
        a set of all genes now in the Mendeliome (and green)
        which were absent in the given panel version
    """

    old: PanelData = {'metadata': [], 'genes': {}}
    get_panel_green(old, old_data={}, version=old_version, forbidden_genes=forbidden)

    return current_genes.difference(set(old['genes'].keys()))


def overwrite_new_status(gene_dict: PanelData, new_genes: set[str]):
    """
    ignores any previous notion of new, replaces it with a manually assigned one

    Args:
        gene_dict ():
        new_genes (): set these genes as new
    """

    panel_id = DEFAULT_PANEL

    for gene, gene_data in gene_dict['genes'].items():
        if gene in new_genes:
            gene_data['new'] = [panel_id]
        else:
            gene_data['new'] = []


@click.command()
@click.option('--panels', help='JSON of per-participant panels')
@click.option('--out_path', required=True, help='destination for results')
@click.option('--dataset', default=None, help='dataset to use, optional')
def main(panels: str | None, out_path: str, dataset: str | None = None):
    """
    if present, reads in any prior reference data
    if present, reads additional panels to use
    queries panelapp for each panel in turn, aggregating results
    optional attribute dataset used for running in a pipeline context

    Args:
        panels (): file containing per-participant panels
        out_path (): where to write the results out to
        dataset (): optional dataset to use
    """

    logging.info('Starting PanelApp Query Stage')

    old_data: dict = {}

    # make responsive to config
    twelve_months = None

    # find and extract this dataset's portion of the config file
    # set the Forbidden genes (defaulting to an empty set)
    forbidden_genes = read_json_from_path(
        get_cohort_config(dataset).get('forbidden', 'missing'), set()
    )

    # historic data overrides default 'previous' list for cohort
    # open to discussing order of precedence here
    if old_file := find_latest_file(start='panel_'):
        logging.info(f'Grabbing legacy panel data from {old_file}')
        old_data: dict = read_json_from_path(old_file, default=old_data)

    elif previous := get_cohort_config(dataset).get('gene_prior'):
        logging.info(f'Reading legacy data from {previous}')
        old_data: dict = read_json_from_path(previous, default=old_data)

    else:
        twelve_months = True

    # are there any genes to skip from the Mendeliome? i.e. only report
    # if in a specifically phenotype-matched panel
    remove_from_core: list[str] = get_config()['panels'].get('require_pheno_match', [])
    logging.info(f'Genes to remove from Mendeliome: {",".join(remove_from_core)!r}')

    # set up the gene dict
    gene_dict: PanelData = {'metadata': [], 'genes': {}}

    # first add the base content
    get_panel_green(
        gene_dict,
        old_data=old_data,
        blacklist=remove_from_core,
        forbidden_genes=forbidden_genes,
    )

    # store the list of genes currently on the core panel
    if twelve_months:
        twelve_months = set(gene_dict['genes'].keys())

    # if participant panels were provided, add each of those to the gene data
    panel_list = set()
    if panels is not None:
        panel_list = read_panels_from_participant_file(panels)
        logging.info(f'Phenotype matched panels: {", ".join(map(str, panel_list))}')

    # now check if there are cohort-wide override panels
    if extra_panels := get_cohort_config(dataset).get('cohort_panels'):
        logging.info(f'Cohort-specific panels: {", ".join(map(str, extra_panels))}')
        panel_list.update(extra_panels)

    for panel in panel_list:

        # skip mendeliome - we already queried for it
        if panel == DEFAULT_PANEL:
            continue

        get_panel_green(
            gene_dict=gene_dict,
            panel_id=panel,
            old_data=old_data,
            forbidden_genes=forbidden_genes,
        )

    # now get the best MOI
    get_best_moi(gene_dict['genes'])

    # if we didn't have prior reference data, scrub down new statuses
    # new_genes can be empty as a result of a successful query
    if twelve_months:

        old_version = find_core_panel_version()
        if old_version is None:
            raise ValueError('Could not find a version from 12 months ago')
        logging.info(
            f'No prior data found, running panel diff vs. panel version {old_version}'
        )
        new_gene_set = get_new_genes(twelve_months, old_version)
        overwrite_new_status(gene_dict, new_gene_set)

    # write the output to long term storage
    write_output_json(output_path=out_path, object_to_write=gene_dict)

    save_new_historic(old_data, prefix='panel_')


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )

    main()
