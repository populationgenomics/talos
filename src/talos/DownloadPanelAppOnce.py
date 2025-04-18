"""
Some deployment environments may not have an internet connection.
To solve this, we now support offline access to PanelApp data by querying once for all content, and referencing as a
local file.
"""

import json

import aiohttp
import asyncio
import re
from argparse import ArgumentParser
from dateutil.parser import parse

from loguru import logger

from talos.config import config_retrieve, ConfigError
from talos.models import (
    PanelShort,
    StolenPanelApp,
    StolenPanelAppGene,
    StolenPanelAppGenePanelDetail,
    StolenPanelAppHpo,
)
from talos.utils import (
    get_json_response,
    read_json_from_path,
)


ENTITY_TYPE_CONSTANT = 'entity_type'
GENE_CONSTANT = 'gene'
HPO_RE = re.compile(r'HP:\d+')
ACTIVITY_CONTENT = {'green list (high evidence)', 'expert review green'}

REALLY_OLD = '1970-01-01'
PANELS_ENDPOINT = 'https://panelapp-aus.org/api/v1/panels'
DEFAULT_PANEL = 137

try:
    PANELS_ENDPOINT = config_retrieve(['GeneratePanelData', 'panelapp'], PANELS_ENDPOINT)
    DEFAULT_PANEL = config_retrieve(['GeneratePanelData', 'default_panel'], DEFAULT_PANEL)
except (ConfigError, KeyError):
    logger.warning('Config environment variable TALOS_CONFIG not set, or keys missing, falling back to Aussie PanelApp')

GREEN_TEMPLATE = f'{PANELS_ENDPOINT}/{{id}}/genes/?confidence_level=3'
ACTIVITY_TEMPLATE = f'{PANELS_ENDPOINT}/{{id}}/activities'


class CustomEncoder(json.JSONEncoder):
    """
    to be used as a JSON encoding class
    - replaces all sets with lists
    - replaces dataclass objects with a dictionary of the same
    """

    def default(self, o):
        """
        takes an arbitrary object, and forms a JSON representation
        where the object doesn't have an easy string representation,
        transform to a valid object: set -> list, class -> dict

        Args:
            o (): python object being JSON encoded
        """

        if isinstance(o, set):
            return list(o)
        return json.JSONEncoder.default(self, o)


def get_panels_and_hpo_terms(endpoint: str = PANELS_ENDPOINT) -> dict[int, StolenPanelAppHpo]:
    """
    query panelapp, collect each panel by its HPO terms
    this dictionary structure is the basis for inserting more information later

    Args:
        endpoint (str): URL for panels

    Returns:
        dict: {panel_ID: StolenPanelAppHpo(**)}
    """

    panels_by_hpo: dict[int, StolenPanelAppHpo] = {}

    while True:
        endpoint_data = get_json_response(endpoint)
        for panel in endpoint_data['results']:
            # create an object to hold this
            panel_hpo = StolenPanelAppHpo(disease_group=panel['disease_group'])
            # can be split over multiple strings, so join then search
            relevant_disorders = ' '.join(panel['relevant_disorders'] or [])
            for match in re.findall(HPO_RE, relevant_disorders):
                panel_hpo.hpo_terms.add(match)

            panels_by_hpo[int(panel['id'])] = panel_hpo

        # cycle through additional pages
        if endpoint := endpoint_data['next']:
            continue
        break

    return panels_by_hpo


def parse_panel_activity(panel_activity: list[dict]) -> dict[str, str]:
    """
    reads in the panel activity dictionary, and for each green entity, finds the date at
    which the entity obtained a Green rating

    Args:
        panel_activity (list[dict]):

    Returns:
        dict, mapping gene symbol to the date it was first graded green (high evidence)
    """

    return_dict: dict[str, str] = {}

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
        creation = parse(activity_entry['created'], ignoretz=True).strftime('%Y-%m-%d')

        # store it
        return_dict[gene_name] = creation

    return return_dict


def parse_panel(
    panel_data: dict,
    panel_activities: list[dict],
    ensg_dict: dict[str, str] | None = None,
    symbol_dict: dict[str, str] | None = None,
) -> dict:
    """

    Args:
        panel_data ():
        panel_activities ():
        ensg_dict (dict): mapping Ensembl IDs to gene symbols, based on MANE data
        symbol_dict (dict): mapping gene symbols to Ensembl IDs, based on MANE data

    Returns:

    """

    # this will contain a range of bits, indexed on ENSG
    panel_gene_content: dict = {}

    green_dates = parse_panel_activity(panel_activities)
    # iterate over the genes in this panel result
    for gene in panel_data['results']:
        if gene['entity_type'] != 'gene':
            continue

        symbol = gene.get('entity_name')

        ensg = None
        mane_ensg = symbol_dict.get(symbol, '') if symbol_dict else ''

        # for some reason the build is capitalised oddly in panelapp, so lower it
        for build, content in gene['gene_data']['ensembl_genes'].items():
            if build.lower() == 'grch38':
                # the ensembl version may alter over time, but will be singular
                ensembl_data = content[next(iter(content.keys()))]
                ensg = ensembl_data['ensembl_id']

        # no ENSG at all, skip completely
        if not (mane_ensg or ensg):
            logger.info(f'Gene {symbol}/{ensg} removed for lack of chrom or ENSG annotation')
            continue

        exact_moi = gene.get('mode_of_inheritance', 'unknown').lower()

        for each_ensg in (ensg, mane_ensg):
            panel_gene_content[each_ensg] = {
                'symbol': symbol,
                'mane_symbol': ensg_dict.get(each_ensg, '') if ensg_dict else '',
                'moi': exact_moi,
                'green_date': green_dates.get(symbol, REALLY_OLD),
            }

    return panel_gene_content


async def get_single_panel(session: aiohttp.ClientSession, panel_id: int) -> dict:
    """

    Args:
        session (aiohttp.ClientSession):
        panel_id ():

    Returns:

    """

    async with session.get(GREEN_TEMPLATE.format(id=panel_id)) as resp:
        reponse = await resp.json()
        return {panel_id: reponse}


async def get_single_panel_activities(session: aiohttp.ClientSession, panel_id: int) -> dict:
    """

    Args:
        session (aiohttp.ClientSession):
        panel_id ():

    Returns:

    """

    async with session.get(ACTIVITY_TEMPLATE.format(id=panel_id)) as resp:
        reponse = await resp.json()
        return {panel_id: reponse}


async def get_all_known_panels(panel_ids: set[int], activities: bool = False) -> dict[int, dict | list]:
    """
    take all the panel IDs, asynchronously query for them
    if panelapp dies it dies

    Args:
        panel_ids ():
        activities (bool): if True, get the activity log for this panel

    Returns:
        the per-panel gene details, or activity log, depending on the activities flag
    """

    tasks = []

    async with aiohttp.ClientSession() as session:
        for panel_id in panel_ids:
            if activities:
                tasks.append(asyncio.ensure_future(get_single_panel_activities(session, panel_id)))
            else:
                tasks.append(asyncio.ensure_future(get_single_panel(session, panel_id)))

        all_panel_details = await asyncio.gather(*tasks)

    return {int(pid): data for panel in all_panel_details for pid, data in panel.items()}


def reorganise_mane_data(mane_path: str) -> tuple[dict[str, str], dict[str, str]]:
    """
    takes the dictionary of MANE data and reorganises into 2 dictionaries
    - one indexed on the gene symbol
    - one on the Ensembl ID

    Args:
        mane_path ():

    Returns:

    """

    raw_mane_data = read_json_from_path(mane_path)
    if not raw_mane_data:
        raise ValueError(f'MANE data not found at {mane_path}')

    ensg_as_primary: dict[str, str] = {}
    symbol_as_primary: dict[str, str] = {}

    for tx_data in raw_mane_data.values():
        symbol = tx_data['symbol']
        ensg = tx_data['ensg']

        ensg_as_primary[ensg] = symbol
        symbol_as_primary[symbol] = ensg

    return ensg_as_primary, symbol_as_primary


def cli_main():
    logger.info('Starting PanelApp parsing')
    parser = ArgumentParser()
    parser.add_argument('--output', help='Where to write Panel data', required=True)
    parser.add_argument('--mane', help='MANE JSON data', default=None)
    args = parser.parse_args()
    main(output=args.output, mane_path=args.mane)


def main(output: str, mane_path: str = None):
    """
    query PanelApp - get EVERYTHING

    Args:
        output (str): path to an output destination
        mane_path (str): path to a MANE JSON file, optional
    """

    # set up a collection object
    collected_panel_data = StolenPanelApp(
        hpos=get_panels_and_hpo_terms(),
    )

    # with open('panelhpo.json', 'w', encoding='utf-8') as handle:
    #     handle.write(collected_panel_data.model_dump_json(indent=4))

    all_panels = set(collected_panel_data.hpos.keys())

    collected_panel_data = read_json_from_path('panelhpo.json', return_model=StolenPanelApp)

    all_panel_data: dict[int, dict] = asyncio.run(get_all_known_panels(all_panels))

    # with open('alldata.json', 'w', encoding='utf-8') as handle:
    #     json.dump(all_panel_data, handle, cls=CustomEncoder)
    #
    # with open('alldata.json', encoding='utf-8') as handle:
    #     all_panel_data = json.load(handle)

    all_panel_activities: dict[int, list] = asyncio.run(get_all_known_panels(all_panels, activities=True))

    # # with open('allacts.json', 'w', encoding='utf-8') as handle:
    # #     json.dump(all_panel_activities, handle, cls=CustomEncoder)
    #
    # with open('allacts.json', encoding='utf-8') as handle:
    #     all_panel_activities = json.load(handle)

    if mane_path:
        ensg_dict, symbol_dict = reorganise_mane_data(mane_path)
    else:
        ensg_dict, symbol_dict = None, None

    # iterate over the gathered panels
    for panel_id, panel_data in all_panel_data.items():
        if not panel_data['results']:
            logger.warning(f'No Green Genes on panel {panel_id}')
            continue

        logger.info(f'Processing panel {panel_id}')

        # get the activity log for this panel
        panel_activities = all_panel_activities[panel_id]

        # parse the data & activities
        parsed_panel_data = parse_panel(
            panel_data=panel_data,
            panel_activities=panel_activities,
            ensg_dict=ensg_dict,
            symbol_dict=symbol_dict,
        )

        one_panel_detail = panel_data['results'][0]['panel']

        collected_panel_data.versions.append(
            PanelShort(
                id=panel_id,
                name=one_panel_detail['name'],
                version=one_panel_detail['version'],
            )
        )

        for gene, gene_data in parsed_panel_data.items():
            # already seen - update some attributes
            if prev_gene_data := collected_panel_data.genes.get(gene):
                prev_gene_data.panels[panel_id] = StolenPanelAppGenePanelDetail(
                    moi=gene_data['moi'],
                    date=gene_data['green_date'],
                )
                # update if previous wasn't populated
                prev_gene_data.mane_symbol = prev_gene_data.mane_symbol or gene_data['mane_symbol']

            else:
                collected_panel_data.genes[gene] = StolenPanelAppGene(
                    symbol=gene_data['symbol'],
                    mane_symbol=gene_data['mane_symbol'],
                    panels={
                        panel_id: StolenPanelAppGenePanelDetail(
                            moi=gene_data['moi'],
                            date=gene_data['green_date'],
                        )
                    },
                )

    with open(output, 'w') as output_file:
        output_file.write(collected_panel_data.model_dump_json(indent=4))


if __name__ == '__main__':
    cli_main()
