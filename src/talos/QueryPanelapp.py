"""
Poll PanelApp for all the Panel contents relevant to this analysis
Start by querying for a 'base' panel (either 137; The Mendeliome, or as specified in the config file)
After the base panel, optionally iterate over all additional panels
Write the data out as a PanelApp object model
"""

from argparse import ArgumentParser
from datetime import datetime
from dateutil.parser import parse
from dateutil.utils import today

from loguru import logger

from talos.config import config_retrieve
from talos.models import PanelApp, PanelDetail, PanelShort, PhenotypeMatchedPanels
from talos.utils import (
    ORDERED_MOIS,
    get_json_response,
    get_simple_moi,
    read_json_from_path,
    parse_mane_json_to_dict,
)


# global variables for PanelApp interaction
PANELAPP_HARD_CODED_DEFAULT = 'https://panelapp-aus.org/api/v1/panels'
# numerical ID of the Mendeliome in PanelApp Australia
PANELAPP_HARD_CODED_BASE_PANEL = 137

try:
    PANELAPP_BASE = config_retrieve(['GeneratePanelData', 'panelapp'], PANELAPP_HARD_CODED_DEFAULT)
    DEFAULT_PANEL = config_retrieve(['GeneratePanelData', 'default_panel'], PANELAPP_HARD_CODED_BASE_PANEL)
except (FileNotFoundError, KeyError):
    logger.warning('Config environment variable TALOS_CONFIG not set, falling back to Aussie PanelApp')
    PANELAPP_BASE = PANELAPP_HARD_CODED_DEFAULT
    DEFAULT_PANEL = PANELAPP_HARD_CODED_BASE_PANEL

ENTITY_TYPE_CONSTANT = 'entity_type'
GENE_CONSTANT = 'gene'
TODAY = today()
REALLY_OLD = parse('1970-01-01')
ACTIVITY_CONTENT = {'green list (high evidence)', 'expert review green'}
DEFAULT_MOI_FOR_MANUAL_GENES = 'Mono_And_Biallelic'
DEFAULT_PANEL_ID = 0


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
    logger.info(f'{panel_name} version: {panel_version}')

    return panel_name, panel_version, panel_genes


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


def get_best_moi(gene_dict: dict):
    """
    From the collected set of all MOIs, take the most lenient
    If set was empty (i.e. no specific MOI found) find the default

    Default is found autosome/sex chrom aware

    Args:
        gene_dict (): the 'genes' index of the collected dict
    """

    for content in gene_dict.values():
        # accept the simplest MOI
        simplified_mois = get_simple_moi(content.all_moi, chrom=content.chrom)

        # force a combined MOI here
        if 'Biallelic' in simplified_mois and 'Monoallelic' in simplified_mois:
            content.moi = 'Mono_And_Biallelic'

        else:
            # take the more lenient of the gene MOI options
            content.moi = sorted(simplified_mois, key=lambda x: ORDERED_MOIS.index(x))[0]


def update_moi_from_config(
    gene_dict: PanelApp,
    add_genes: list[dict[str, str]],
    gene_json: str | None = None,
) -> PanelApp:
    """
    Add additional genes from config
    Update the MOI for a gene or genes based on a dictionary of gene: moi

    Args:
        gene_dict (PanelApp): the gene data to update
        add_genes (list[dict]): entities to add
        gene_json (str | None): path to the symbol: ensg dictionary in JSON form
    """

    logger.info(f'Overriding MOI for specific genes: {add_genes}')

    # read the MANE file into a dictionary
    mane_lookup = parse_mane_json_to_dict(gene_json) if gene_json else {}

    for gene_data in add_genes:
        # this field is mandatory
        ensg = gene_data['ensg']

        # not essential, but if we have a symbol lookup, use it
        # if it wasn't provided, use the lookup, if it's not in the lookup, use ENSG
        symbol = gene_data.get('symbol', mane_lookup.get(ensg, ensg))

        # this field is optional
        moi = gene_data.get('moi', DEFAULT_MOI_FOR_MANUAL_GENES)

        if moi not in ORDERED_MOIS:
            raise ValueError(f'{moi} for {ensg} is not a valid MOI, choose from {", ".join(ORDERED_MOIS)}')

        # if this is a gene we already know about, try and update the MOI
        # and add the custom panel ID to associated panels
        if ensg in gene_dict.genes:
            current_moi = gene_dict.genes[ensg].moi
            if current_moi == moi:
                logger.info(f'{ensg} MOI was already {moi}, no change made')
            else:
                logger.info(f'{ensg} MOI was updated from {current_moi} to {moi}')
                gene_dict.genes[ensg].moi = moi
            gene_dict.genes[ensg].panels.add(DEFAULT_PANEL_ID)

        else:
            # if this is a new gene, add it to the dictionary
            gene_dict.genes[ensg] = PanelDetail(
                symbol=symbol,
                all_moi={moi},
                moi=moi,
                panels={DEFAULT_PANEL_ID},
                chrom=gene_data.get('chrom', 'chrUnknown'),
            )
    return gene_dict


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--input', help='JSON of per-participant panels', default=None)
    parser.add_argument('--output', required=True, help='destination for results')
    parser.add_argument('--gene_json', help='path to gene symbol JSON', default=None)
    args = parser.parse_args()
    main(panels=args.input, out_path=args.output, gene_json_path=args.gene_json)


def main(panels: str | None, out_path: str, gene_json_path: str | None = None):
    """
    Queries PanelApp for all the gene panels to use in the current analysis
    queries panelapp for each panel in turn, aggregating results

    Args:
        panels (): file containing per-participant panels
        out_path (): where to write the results out to
        gene_json_path (): path to gene symbol JSON, used to obtain a lookup of symbol to ENSG ID
    """

    logger.info('Starting PanelApp Query Stage')

    # set the Forbidden genes (defaulting to an empty set)
    forbidden_genes = config_retrieve(['GeneratePanelData', 'forbidden_genes'], set())

    # are there any genes to skip from the Mendeliome? i.e. only report if in a specifically phenotype-matched panel
    remove_from_core: list[str] = config_retrieve(['GeneratePanelData', 'require_pheno_match'], [])
    logger.info(f'Genes to remove from Mendeliome: {",".join(remove_from_core)!r}')

    # set up the gene dict
    gene_dict = PanelApp(genes={})

    # first add the base content
    logger.info('Getting Base Panel')
    get_panel(gene_dict, blacklist=remove_from_core, forbidden_genes=forbidden_genes)

    # if participant panels were provided, add each of those to the gene data
    panel_list: set[int] = set()
    if panels is not None:
        logger.info('Reading participant panels')
        hpo_panel_object = read_json_from_path(panels, return_model=PhenotypeMatchedPanels)
        panel_list = hpo_panel_object.all_panels
        logger.info(f'Phenotype matched panels: {", ".join(map(str, panel_list))}')

    # now check if there are cohort-wide override panels
    if extra_panels := config_retrieve(['GeneratePanelData', 'forced_panels'], False):
        logger.info(f'Cohort-specific panels: {", ".join(map(str, extra_panels))}')
        panel_list.update(extra_panels)

    for panel in panel_list:
        # skip mendeliome - we already queried for it
        if panel == DEFAULT_PANEL:
            continue
        logger.info(f'Getting Panel {panel}')
        get_panel(gene_dict=gene_dict, panel_id=panel, forbidden_genes=forbidden_genes)

    # now get the best MOI, and update the entities in place
    get_best_moi(gene_dict.genes)

    # try and pull from config, default to an empty list
    if gene_data := config_retrieve(['GeneratePanelData', 'manual_overrides', 'genes'], []):
        gene_dict = update_moi_from_config(gene_dict, gene_data, gene_json_path)

    # write the output to long term storage
    with open(out_path, 'w') as out_file:
        out_file.write(PanelApp.model_validate(gene_dict).model_dump_json(indent=4))


if __name__ == '__main__':
    cli_main()
