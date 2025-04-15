"""
Some deployment environments may not have an internet connection.
To solve this, we now support offline access to PanelApp data by querying once for all content, and referencing as a
local file.
"""

import re
from argparse import ArgumentParser
from collections import defaultdict

from datetime import datetime
from dateutil.parser import parse
from dateutil.utils import today

from loguru import logger

from talos.config import config_retrieve, ConfigError
from talos.models import PanelApp, PanelDetail, PanelShort, PhenotypeMatchedPanels
from talos.utils import (
    get_json_response,
    read_json_from_path,
    parse_mane_json_to_dict,
)


# # TODO
# https://panelapp-aus.org/api/v1/panels/
#
# then for each panel:
#
# https://panelapp-aus.org/api/v1/panels/76/genes/?confidence_level=3
#
# - save every panel Id and its HPO terms
# - save every unique gene and the mapping of panel ID: MOI


HPO_RE = re.compile(r'HP:\d+')

PANELS_ENDPOINT = 'https://panelapp-aus.org/api/v1/panels'
DEFAULT_PANEL = 137

try:
    PANELS_ENDPOINT = config_retrieve(['GeneratePanelData', 'panelapp'], PANELS_ENDPOINT)
    DEFAULT_PANEL = config_retrieve(['GeneratePanelData', 'default_panel'], DEFAULT_PANEL)
except (ConfigError, KeyError):
    logger.warning('Config environment variable TALOS_CONFIG not set, or keys missing, falling back to Aussie PanelApp')


def get_panels_and_hpo_terms(endpoint: str = PANELS_ENDPOINT) -> dict[int, dict[str, str | set[str]]]:
    """
    query panelapp, collect each panel by its HPO terms
    this dictionary structure is the basis for inserting more information later

    Args:
        endpoint (str): URL for panels

    Returns:
        dict: {panel_ID: {'hpo_terms': [HPO_Term, HPO_Term], }}
    """

    panels_by_hpo: dict[int, dict] = defaultdict(dict)

    while True:
        endpoint_data = get_json_response(endpoint)
        for panel in endpoint_data['results']:
            panels_by_hpo[int(panel['id'])] = {
                'hpo_terms': set(),
                'disease_group': panel['disease_group'],
                'name': panel['name'],
                'version': panel['version'],
            }
            # can be split over multiple strings
            relevant_disorders = ' '.join(panel['relevant_disorders'] or [])
            for match in re.findall(HPO_RE, relevant_disorders):
                panels_by_hpo[int(panel['id'])]['hpo_terms'].add(match)

        # cycle through additional pages
        if endpoint := endpoint_data['next']:
            continue
        break

    return dict(panels_by_hpo)


def parse_panel_activity(panel_activity: list[dict]) -> dict[str, datetime]:
    """
    reads in the panel activity dictionary, and for each green entity, finds the date at
    which the entity obtained a Green rating

    Args:
        panel_activity (list[dict]):

    Returns:
        dict, mapping gene symbol to the date it was first graded green (high evidence)
    """

    return_dict: dict[str, datetime] = {}

    # do some stuff
    for activity_entry in panel_activity:
        # only interested in genes at the moment
        if activity_entry.get(ENTITY_TYPE_CONSTANT) != GENE_CONSTANT:
            continue

        # get the name of the gene
        gene_name = activity_entry['entity_name']

        # check for relevant text
        lower_text = activity_entry['text'].lower()
        if not any(each_string in lower_text for each_string in ACTIVITY_CONTENT):
            continue

        # find the event date for this activity entry
        creation = parse(activity_entry['created'], ignoretz=True)

        # store it
        return_dict[gene_name] = creation
    return return_dict


def get_panel(
    gene_dict: PanelApp,
    panel_id: int = DEFAULT_PANEL,
    blacklist: list[str] | None = None,
    forbidden_genes: set[str] | None = None,
):
    """
    Takes a panel number, and pulls all GRCh38 gene details from PanelApp
    For each gene, keep the MOI, symbol, ENSG (where present)

    Args:
        gene_dict (): PanelApp obj to continue populating
        panel_id (): specific panel or 'base' (e.g. 137)
        blacklist (): list of symbols/ENSG IDs to remove from this panel
        forbidden_genes (set[str]): genes to remove for this cohort
    """

    if blacklist is None:
        blacklist = []

    if forbidden_genes is None:
        forbidden_genes = set()

    panel_name, panel_version, panel_genes = request_panel_data(f'{PANELAPP_BASE}/{panel_id}/')

    # get the activity log for this panel
    panel_activity = get_json_response(f'{PANELAPP_BASE}/{panel_id}/activities/')

    green_dates = parse_panel_activity(panel_activity)

    # find the threshold for when a gene should be treated as recent - new if added within this many months
    # by default we're falling back to 6 months, just so we don't fail is this is absent in config
    recent_months = config_retrieve(['GeneratePanelData', 'within_x_months'], 6)

    # add metadata for this panel & version
    gene_dict.metadata.append(PanelShort(name=panel_name, version=panel_version, id=panel_id))

    # iterate over the genes in this panel result
    for gene in panel_genes:
        symbol = gene.get('entity_name')

        # how long ago was this added to this panel? If within the last X months, treat as new
        # if we didn't find an acceptable date from the API, fall back on REALLY_OLD (never recent)
        # relativedelta is complete ass for this test, rewriting manually here
        # for posterity, relativedelta in dateutil does this calculation, then overwrites it with a
        # non year-aware version, which is a bit of a mess IMO
        # the dateutil result between March 2023 and September 2024 is 6 months, which is incorrect
        # for this purpose as it ignores the 12 full months between the two dates
        event_datetime = green_dates.get(symbol, REALLY_OLD)
        months = (TODAY.year - event_datetime.year) * 12 + (TODAY.month - event_datetime.month)
        new_gene = months < recent_months

        # only retain green genes
        if gene['confidence_level'] != '3' or gene['entity_type'] != 'gene' or symbol in forbidden_genes:
            continue

        ensg = None
        chrom = None

        # TODO stop picking up ENSG from here
        # for some reason the build is capitalised oddly in panelapp
        # at least one entry doesn't have an ENSG annotation
        for build, content in gene['gene_data']['ensembl_genes'].items():
            if build.lower() == 'grch38':
                # the ensembl version may alter over time, but will be singular
                ensembl_data = content[next(iter(content.keys()))]
                ensg = ensembl_data['ensembl_id']
                chrom = ensembl_data['location'].split(':')[0]

        if chrom is None:
            logger.info(f'Gene {symbol}/{ensg} removed from {panel_name} for lack of chrom annotation')
            continue

        if ensg is None or ensg in blacklist or symbol in blacklist or ensg in forbidden_genes:
            logger.info(f'Gene {symbol}/{ensg} removed from {panel_name}')
            continue

        exact_moi = gene.get('mode_of_inheritance', 'unknown').lower()

        # either update or add a new entry
        if ensg in gene_dict.genes:
            this_gene = gene_dict.genes[ensg]

            # now we find it on this panel
            this_gene.panels.add(panel_id)

            # add this moi to the set
            this_gene.all_moi.add(exact_moi)

            # if this is a recent addition to the panel, add this panel to the 'new' attr
            if new_gene:
                this_gene.new.add(panel_id)

        else:
            # save the entity into the final dictionary
            gene_dict.genes[ensg] = PanelDetail(
                symbol=symbol,
                all_moi={exact_moi},
                new={panel_id} if new_gene else set(),
                panels={panel_id},
                chrom=chrom,
            )


def cli_main():
    logger.info('Starting HPO~Panel matching')
    parser = ArgumentParser()
    parser.add_argument('--output', help='Where to write Panel data', required=True)
    args = parser.parse_args()
    main(output=args.output)


def main(output: str):
    """
    query PanelApp - get EVERYTHING

    Args:
        output (str): path to an output destination
    """
    print(output)
    panels_and_hpo_terms = get_panels_and_hpo_terms()


if __name__ == '__main__':
    cli_main()
