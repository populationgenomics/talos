"""
Combines the two behaviours of
- GeneratePanelData: use HPO terms attached to participants to match panels to families
- QueryPanelApp: Use the PanelApp API to get all genes and panels relevant to the analysis

Takes as input:
- The downloaded PanelApp data
- Optionally a HPO obo file & participant phenopackets. These last two should be coupled - both either absent or present
"""

from argparse import ArgumentParser
from collections import defaultdict
from string import punctuation

import pendulum
import phenopackets.schema.v2 as pps2
from google.protobuf.json_format import ParseDict
from loguru import logger
from networkx import dfs_successors
from networkx.exception import NetworkXError
from obonet import read_obo

from talos.config import config_retrieve
from talos.models import (
    PanelApp,
    PanelShort,
    DownloadedPanelApp,
    PanelDetail,
    ParticipantHPOPanels,
    PhenoPacketHpo,
)
from talos.utils import read_json_from_path

PANELAPP_BASE_PANEL = 137
X_CHROMOSOME = {'X'}
# most lenient to most conservative
# usage = if we have two MOIs for the same gene, take the broadest
ORDERED_MOIS = ['Mono_And_Biallelic', 'Monoallelic', 'Hemi_Mono_In_Female', 'Hemi_Bi_In_Female', 'Biallelic']
IRRELEVANT_MOI = {'unknown', 'other'}

# we consider a gene new, and worth flagging, if it was made Green in a panel within this time frame
WITHIN_X_MONTHS = 6

try:
    DEFAULT_PANEL = config_retrieve(['GeneratePanelData', 'default_panel'], PANELAPP_BASE_PANEL)
    WITHIN_X_MONTHS = config_retrieve(['GeneratePanelData', 'within_x_months'], WITHIN_X_MONTHS)
except KeyError:
    logger.warning('Config environment variable TALOS_CONFIG not set, falling back to Aussie PanelApp')
    DEFAULT_PANEL = PANELAPP_BASE_PANEL

# create a datetime threshold
NEW_THRESHOLD = pendulum.now().subtract(months=WITHIN_X_MONTHS)

MOI_FOR_CUSTOM_GENES = 'Mono_And_Biallelic'
CUSTOM_PANEL_ID = 0


def cli_main():
    logger.info('Starting HPO~Panel matching')
    parser = ArgumentParser()
    parser.add_argument('--input', help='Pre-Downloaded PanelApp data', required=True)
    parser.add_argument('--output', help='Path to write JSON output to', required=True)
    parser.add_argument('--cohort', help='GA4GH Cohort/PhenoPacket File', required=True)
    parser.add_argument('--hpo', help='Localised copy of HPO obo file', required=False)
    args = parser.parse_args()
    main(panel_data=args.input, output_file=args.output, cohort_file=args.cohort, hpo_file=args.hpo)


def extract_participant_data_from_cohort(cohort: pps2.Cohort) -> tuple[PanelApp, set[str]]:
    """
    read the extended pedigree file, pull out family details and HPO terms

    Args:
        cohort (str): GA4GH Cohort/PhenoPacket file

    Returns:
        set of all HPO terms
    """

    # generate an empty result holder
    shell = PanelApp()

    all_hpos: set[str] = set()

    for member in cohort.members:
        shell.participants[member.id] = ParticipantHPOPanels(
            external_id=member.subject.alternate_ids[0] if member.subject.alternate_ids else member.id,
            family_id=member.subject.id,
            hpo_terms=[PhenoPacketHpo(id=hpo.type.id, label=hpo.type.label) for hpo in member.phenotypic_features],
            panels={DEFAULT_PANEL},
        )
        all_hpos.update({hpo.type.id for hpo in member.phenotypic_features})

    return shell, all_hpos


def match_hpos_to_panels(
    hpo_panel_map: dict[int, list[PhenoPacketHpo]],
    hpo_file: str,
    all_hpos: set[str],
) -> dict[str, set[int]]:
    """
    take the HPO terms from the participant metadata, and match to panels

    Args:
        hpo_panel_map (dict): dict of panel IDs to hpo details
        hpo_file (str): path to an obo file containing HPO tree
        all_hpos (set[str]): collection of all unique HPO terms

    Returns:
        a dictionary linking all HPO terms to a corresponding set of Panel IDs
    """

    hpo_graph = read_obo(hpo_file, ignore_obsolete=False)

    # re-index the HPO terms to the panel IDs
    panel_per_hpo = defaultdict(set)
    for panel_id, hpo_terms in hpo_panel_map.items():
        for phenopacket_hpo in hpo_terms:
            panel_per_hpo[phenopacket_hpo.id].add(panel_id)

    # create a fresh object to hold all matched hpos to panels
    hpo_to_panels = defaultdict(set)

    # cycle through all HPO terms in the cohort
    # for each term chase it back to the HPO ontology root
    # if there are any hits between the HPO chain and a panel HPO, retain that connection
    for hpo_string in all_hpos:
        # identify all HPO terms back to the ontology root
        try:
            successor_hpo_terms = set(dfs_successors(hpo_graph, hpo_string))

            for hpo_term in successor_hpo_terms:
                if hpo_term in panel_per_hpo:
                    hpo_to_panels[hpo_string].update(panel_per_hpo[hpo_term])
        except (KeyError, NetworkXError):
            # if the HPO term is not in the graph, skip it
            logger.warning(f'HPO term {hpo_string} not found in HPO graph')

    return dict(hpo_to_panels)


def match_participants_to_panels(panelapp_data: PanelApp, hpo_panels: dict, cached_panelapp: DownloadedPanelApp):
    """
    take the two maps of Participants: HPOs, and HPO: Panels
    blend the two to find panels per participant

    For each participant, find any HPO terms which were matched to panels
    for each matched term, add the panel(s) to the participant's private set

    Args:
        panelapp_data (PanelApp): CPG ID to phenotype details
        hpo_panels (dict): lookup of panels per HPO term
        cached_panelapp (DownloadedPanelApp): pre-downloaded PanelApp data, steal versions from here
    """

    forced_panels: set[int] = config_retrieve(['GeneratePanelData', 'forced_panels'], set())

    all_panels_in_this_analysis: set[int] = {DEFAULT_PANEL} | forced_panels

    for party_data in panelapp_data.participants.values():
        for hpo_term in party_data.hpo_terms:
            if panel_list := hpo_panels.get(hpo_term.id):
                # add relevant panels for this participant
                party_data.panels.update(panel_list)
                # and add to the collection of all panels
                all_panels_in_this_analysis.update(panel_list)

        if forced_panels:
            # add the forced panels to the participant
            party_data.panels.update(forced_panels)

    for panel in cached_panelapp.versions:
        # check if the panel is in the analysis (i.e. matched to a family)
        if panel.id in all_panels_in_this_analysis:
            # add to the list of panels for this analysis
            panelapp_data.metadata[panel.id] = panel


def get_simple_moi(input_mois: set[str], chrom: str) -> str:
    """
    takes the various PanelApp MOIs, and reduces to a range of cases which can be easily implemented in RD analysis

    Step 1 is finding a simplified MOI for each known PanelApp String
    Step 2 is using the collection of simplified terms to boil this down to one consensus MOI for the gene

    Args:
        input_mois (set[str]): all the MOIs for this gene
        chrom ():
    """

    default = 'Hemi_Bi_In_Female' if chrom in X_CHROMOSOME else 'Biallelic'

    simplified_mois: set[str] = set()

    for input_moi in input_mois:
        # skip over ignore-able MOIs
        if input_moi in IRRELEVANT_MOI:
            continue

        # split each PanelApp MOI into a list of strings
        input_list = input_moi.translate(str.maketrans('', '', punctuation)).split()

        # run a match: case to classify it
        match input_list:
            case ['biallelic', *_additional]:
                simplified_mois.add('Biallelic')
            case ['both', *_additional]:
                simplified_mois.add('Mono_And_Biallelic')
            case ['monoallelic', *_additional]:
                if chrom in X_CHROMOSOME:
                    simplified_mois.add('Hemi_Mono_In_Female')
                else:
                    simplified_mois.add('Monoallelic')
            case ['xlinked', *additional] if 'biallelic' in additional:
                simplified_mois.add('Hemi_Bi_In_Female')
            case ['xlinked', *_additional]:
                simplified_mois.add('Hemi_Mono_In_Female')
            case _:
                continue

    # adda default - solves the all-irrelevant or empty-input cases
    if not simplified_mois:
        return default

    # force a combined MOI here
    if 'Biallelic' in simplified_mois and 'Monoallelic' in simplified_mois:
        return 'Mono_And_Biallelic'

    if all(x in simplified_mois for x in ['Hemi_Bi_In_Female', 'Hemi_Mono_In_Female']):
        return 'Hemi_Mono_In_Female'

    # take the more lenient of the gene MOI options
    return sorted(simplified_mois, key=lambda x: ORDERED_MOIS.index(x))[0]


def fetch_genes_for_panels(panelapp_data: PanelApp, cached_panelapp: DownloadedPanelApp):
    """
    now that we know which panels will be in the analysis, get the corresponding genes and consensus MOI for each

    Args:
        panelapp_data ():
        cached_panelapp ():

    Returns:

    """

    full_set_of_panels: set[int] = set()

    # collect all the panels in the analysis
    for participant_details in panelapp_data.participants.values():
        full_set_of_panels.update(participant_details.panels)

    # now iterate over the genes we know about, and pull them into the panel details object
    for gene_data in cached_panelapp.genes.values():
        # check if this gene is in any of the panels we are interested in - retain the overlap
        if not (panel_intersection := set(gene_data.panels.keys()).intersection(full_set_of_panels)):
            continue

        # collect all the relevant MOIs based on the panels we're using
        gene_mois = {gene_data.panels[panel_id].moi for panel_id in panel_intersection}

        # get the consensus MOI for this gene
        moi = get_simple_moi(gene_mois, gene_data.chrom)

        # find any panels where this gene is new
        new_panels = {
            panel_id
            for panel_id in panel_intersection
            if pendulum.from_format(gene_data.panels[panel_id].date, 'YYYY-MM-DD') > NEW_THRESHOLD
        }

        # add the gene to the panel details object
        panelapp_data.genes[gene_data.ensg] = PanelDetail(
            symbol=gene_data.symbol,
            chrom=gene_data.chrom,
            moi=moi,
            new=new_panels,
            panels=panel_intersection,
        )


def update_moi_from_config(
    panelapp_data: PanelApp,
    add_genes: list[dict[str, str]],
):
    """
    Add additional genes from config
    Update the MOI for a gene or genes based on a dictionary of gene: moi

    Args:
        panelapp_data (PanelApp): the gene data to update
        add_genes (list[dict]): entities to add
    """

    logger.info(f'Updating with custom gene content: {add_genes}')

    panelapp_data.metadata[0] = PanelShort(id=0, name='Custom data from Config')

    for gene_data in add_genes:
        # this field is mandatory
        ensg = gene_data['ensg']

        # this field is optional, falling back to a lenient default
        moi = gene_data.get('moi', MOI_FOR_CUSTOM_GENES)

        if moi not in ORDERED_MOIS:
            raise ValueError(f'{moi} for {ensg} is not a valid MOI, choose from {", ".join(ORDERED_MOIS)}')

        # if this is a gene we already know about, update the MOI and add the custom panel ID to associated panels
        if ensg in panelapp_data.genes:
            logger.info(f'From custom data: setting {ensg} MOI to {moi}')
            panelapp_data.genes[ensg].moi = moi
            panelapp_data.genes[ensg].panels.add(CUSTOM_PANEL_ID)

        else:
            # if this is a new gene, add it to the dictionary
            # if we have a symbol use it, if it wasn't provided, use ENSG
            panelapp_data.genes[ensg] = PanelDetail(
                symbol=gene_data.get('symbol', f'Custom: {ensg}'),
                moi=moi,
                panels={CUSTOM_PANEL_ID},
                chrom=gene_data.get('chrom', 'chrUnknown'),
            )


def remove_blacklisted_genes(panelapp_data: PanelApp):
    """
    remove any genes which are blacklisted in the config
    Args:
        panelapp_data ():
    """

    # if any genes are blacklisted in config, remove them here
    if forbidden_genes := config_retrieve(['GeneratePanelData', 'forbidden_genes'], set()):
        genes_to_remove = set(forbidden_genes).intersection(set(panelapp_data.genes.keys()))
        for gene in genes_to_remove:
            del panelapp_data.genes[gene]


def main(panel_data: str, output_file: str, cohort_file: str | None = None, hpo_file: str | None = None):
    """
    Loads the pre-downloaded PanelApp content

    If a Cohort and HPO obo file are present, we attempt matching of panels to phenotypes

    Based on all the panels in the analysis, we then fetch the relevant genes and their MOIs

    Args:
        panel_data ():
        output_file ():
        cohort_file ():
        hpo_file ():
    """

    cached_panelapp: DownloadedPanelApp = read_json_from_path(panel_data, return_model=DownloadedPanelApp)

    # catch all the HPO terms in the cohort, and grab the participant metadata
    cohort = ParseDict(read_json_from_path(cohort_file), pps2.Cohort())

    # extract participant metadata from the Cohort, and collect each unique HPO term
    panelapp_data, all_hpos = extract_participant_data_from_cohort(cohort=cohort)

    # chuck in the default Mendeliome metadata
    panelapp_data.metadata = {DEFAULT_PANEL: cached_panelapp.versions[DEFAULT_PANEL]}

    if hpo_file:
        logger.info('HPO file present, running panel matching algorithm')

        # match HPO terms to panel IDs
        # returns a lookup of each HPO term in the cohort, and panels it is associated with
        hpo_to_panels = match_hpos_to_panels(
            hpo_panel_map=cached_panelapp.hpos,
            all_hpos=all_hpos,
            hpo_file=hpo_file,
        )

        # associate those panels with participants in the cohort
        match_participants_to_panels(panelapp_data, hpo_to_panels, cached_panelapp)

    # now that we have the panels to use, go get them, and assign a single MOI to each gene
    fetch_genes_for_panels(panelapp_data=panelapp_data, cached_panelapp=cached_panelapp)

    # optionally shove in some extra gene content from configuration as a custom panel
    if custom_content := config_retrieve(['PanelApp', 'manual_overrides'], []):
        update_moi_from_config(panelapp_data, custom_content)

    # if any genes are blacklisted in config, remove them here
    remove_blacklisted_genes(panelapp_data)

    # validate and write using pydantic
    valid_cohort_details = PanelApp.model_validate(panelapp_data)
    with open(output_file, 'w', encoding='utf-8') as handle:
        handle.write(valid_cohort_details.model_dump_json(indent=4))


if __name__ == '__main__':
    cli_main()
