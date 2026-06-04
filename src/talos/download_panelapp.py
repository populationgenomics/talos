"""
Prior to the creation of this script, PanelApp was queried on two separate occasions per run:

- once to get all the panels in PanelApp, and the associated HPO terms for each
- once to get the gene data for each panel

There's a couple of ways in which this is wasteful:
- PanelApp is only updated once per month, but this querying goes on in every run
- the original API query implementation was in series, not parallel, due to server instability, so it was slow

In addition to the above, the core assumption here is that deployment environments will have an internet connection.
That is typically not true of HPC environments, making this a tough workflow to run in those environments.

This script queries for the whole of PanelApp in one go, using a barrage of highly parallelised queries:
- use the panels endpoint to get all the panels, and their associated HPO terms
- use the panel-genes endpoint to pull all the gene content for each panel
- use the panel-activities endpoint to pull the activity log for each panel

Save ALL this data:
- collect a list of all current panel IDs/names/versions
- for each panel, collect the associated phenotypic terms
- for each gene, collect the IDs of the panels it appears on, and the date it was first graded green

Optionally takes a MANE JSON file, which is used to map Ensembl IDs to gene symbols, & vice versa. If supplied:
- attempt to find alternative gene symbols for the ENSG ID
- attempt to find alternative Ensembl IDs for the gene symbol
- record all variations
"""

import asyncio
import re
from argparse import ArgumentParser

import aiohttp
from dateutil.parser import parse
from loguru import logger

from talos.config import ConfigError, config_retrieve
from talos.models import (
    DownloadedPanelApp,
    DownloadedPanelAppGene,
    DownloadedPanelAppGenePanelDetail,
    HpoTerm,
    PanelShort,
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
    DEFAULT_PANEL = config_retrieve(['GeneratePanelData', 'default_panel'], DEFAULT_PANEL)
    PANELS_ENDPOINT = config_retrieve(['GeneratePanelData', 'panelapp'], PANELS_ENDPOINT)
except (ConfigError, KeyError):
    logger.warning('Config environment variable TALOS_CONFIG not set, or keys missing, falling back to Aussie PanelApp')

# if this is a massive result, it returns over a number of pages
PANEL_TEMPLATE_URL = f'{PANELS_ENDPOINT}/{{id}}'
ACTIVITY_TEMPLATE = f'{PANELS_ENDPOINT}/{{id}}/activities'
MITO_BAD = 'MT'
MITO_GOOD = 'M'


def get_panels_and_hpo_terms(endpoint: str = PANELS_ENDPOINT) -> dict[int, list[HpoTerm]]:
    """
    query panelapp, collect each panel by its HPO terms

    Args:
        endpoint (str): URL for panels

    Returns:
        dict: {panel_ID: [HPO term, HPO term]}
    """

    panels_by_hpo: dict[int, list[HpoTerm]] = {}

    while True:
        endpoint_data = get_json_response(endpoint)
        for panel in endpoint_data['results']:
            panel_id = int(panel['id'])

            panels_by_hpo[panel_id] = []

            # can be split over multiple strings, so join then search
            relevant_disorders = ' '.join(panel['relevant_disorders'] or [])
            for match in re.findall(HPO_RE, relevant_disorders):
                panels_by_hpo[panel_id].append(HpoTerm(id=match, label=''))

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
    panel_data: dict[str, str | list[dict]],
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
    """

    # this will contain a range of bits, indexed on ENSG
    panel_gene_content: dict = {}

    green_dates: dict[str, str] = parse_panel_activity(panel_activities)

    # iterate over the genes in this panel result
    for gene in panel_data['genes']:
        # please the linter
        if not isinstance(gene, dict):
            raise TypeError(f'Gene {gene} is not a dict')

        symbol: str = gene['symbol']
        ensg: str = gene['ensg']
        mane_ensg = symbol_dict.get(symbol, '') if symbol_dict else ''

        # no ENSG at all, skip completely
        if not (mane_ensg or ensg):
            logger.info(f'Gene {symbol}/{ensg} removed for lack of chrom or ENSG annotation')
            continue

        for each_ensg in [ensg, mane_ensg]:
            if not each_ensg:
                continue

            panel_gene_content[each_ensg] = {
                'symbol': symbol,
                'chrom': gene['chrom'],
                'mane_symbol': ensg_dict.get(each_ensg, '') if ensg_dict else '',
                'moi': gene['moi'],
                'green_date': green_dates.get(symbol, REALLY_OLD),
                'confidence_level': gene['confidence_level'],
            }

    return panel_gene_content


async def get_single_panel(session: aiohttp.ClientSession, panel_id: int) -> dict[int, dict[str, str | list[dict]]]:
    """
    Async method to return data from a single panel.
    Does most of the initial parsing of panel data to reduce memory footprint.

    Args:
        session: aiohttp ClientSession
        panel_id: int, panel ID to search for

    Returns:
        dict, indexed by panel ID, containing panel genes, name, and version
    """
    panel_url = PANEL_TEMPLATE_URL.format(id=panel_id)
    gene_results: list[dict] = []

    async with session.get(panel_url) as resp:
        response = await resp.json()

        # thin out the results, what do we need?
        panel_name = response['name']
        panel_version = response['version']
        for gene in response['genes']:
            # genes only here for now
            if gene['entity_type'] != 'gene':
                continue

            chrom: str = ''
            ensg: str | None = None

            # for some reason the build is capitalised oddly in panelapp, so lower it
            for build, content in gene['gene_data']['ensembl_genes'].items():
                if build.lower() == 'grch38':
                    # the ensembl version may alter over time, but will be singular
                    ensembl_data = content[next(iter(content.keys()))]
                    ensg = ensembl_data['ensembl_id']
                    chrom = ensembl_data['location'].split(':')[0]

                    # step this down to M, for Hail
                    if chrom == MITO_BAD:
                        chrom = MITO_GOOD

            gene_results.append(
                {
                    'symbol': gene['entity_name'],
                    'chrom': chrom,
                    'ensg': ensg,
                    'moi': gene.get('mode_of_inheritance', 'unknown').lower(),
                    'confidence_level': int(gene['confidence_level']),
                }
            )

    return {panel_id: {'name': panel_name, 'version': panel_version, 'genes': gene_results}}


async def get_single_panel_activities(session: aiohttp.ClientSession, panel_id: int) -> dict:
    """Async method to get activities from a single panel"""

    async with session.get(ACTIVITY_TEMPLATE.format(id=panel_id)) as resp:
        reponse = await resp.json()
        return {panel_id: reponse}


async def get_all_known_panels(panel_ids: set[int], activities: bool = False) -> dict:
    """Take all the panel IDs, asynchronously query for them. If panelapp dies it dies."""

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


def main(output: str, mane_path: str | None = None):
    """
    query PanelApp - get EVERYTHING

    Args:
        output (str): path to an output destination
        mane_path (str): path to a MANE JSON file, optional
    """

    # set up a collection object - loaded method execution
    collected_panel_data = DownloadedPanelApp(hpos=get_panels_and_hpo_terms())

    all_panels = set(collected_panel_data.hpos.keys())

    async def _fetch_all() -> tuple[dict, dict]:
        return await asyncio.gather(
            get_all_known_panels(all_panels),
            get_all_known_panels(all_panels, activities=True),
        )

    all_panel_data, all_panel_activities = asyncio.run(_fetch_all())

    if mane_path:
        ensg_dict, symbol_dict = reorganise_mane_data(mane_path)
    else:
        ensg_dict, symbol_dict = None, None

    zero_green_panels: list[int] = []

    # iterate over the gathered panels
    for panel_id, panel_data in all_panel_data.items():
        if not panel_data['genes']:
            logger.warning(f'No Genes on panel {panel_id}')
            zero_green_panels.append(panel_id)
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

        collected_panel_data.versions.append(
            PanelShort(
                id=panel_id,
                name=panel_data['name'],
                version=panel_data['version'],
            ),
        )

        for gene, gene_data in parsed_panel_data.items():
            # already seen - update some attributes
            if prev_gene_data := collected_panel_data.genes.get(gene):
                prev_gene_data.panels[panel_id] = DownloadedPanelAppGenePanelDetail(
                    moi=gene_data['moi'],
                    date=gene_data['green_date'],
                    confidence=gene_data['confidence_level'],
                )
                # update if previous wasn't populated
                prev_gene_data.mane_symbol = prev_gene_data.mane_symbol or gene_data['mane_symbol']

            else:
                collected_panel_data.genes[gene] = DownloadedPanelAppGene(
                    chrom=gene_data['chrom'],
                    symbol=gene_data['symbol'],
                    ensg=gene,
                    mane_symbol=gene_data['mane_symbol'],
                    panels={
                        panel_id: DownloadedPanelAppGenePanelDetail(
                            moi=gene_data['moi'],
                            date=gene_data['green_date'],
                            confidence=gene_data['confidence_level'],
                        ),
                    },
                )

    # strip out any panels with no green genes on, so they're not considered for HPO matches
    for panel_id in zero_green_panels:
        logger.info(f'Removing panel {panel_id} from hpo matching - no green genes')
        del collected_panel_data.hpos[panel_id]

    with open(output, 'w') as output_file:
        output_file.write(collected_panel_data.model_dump_json(indent=4))


if __name__ == '__main__':
    cli_main()
