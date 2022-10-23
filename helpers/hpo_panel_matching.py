"""
all methods relating to:
- querying metamist for participants and their HPO terms
- querying PanelApp for all panels with HPO terms
- loose-matching all HPO terms to relevant panels
- collecting a list of relevant panels for all participants
"""


import logging
import json
import sys
import re
from argparse import ArgumentParser
from collections import defaultdict

import networkx
from obonet import read_obo
from sample_metadata.apis import SeqrApi

from cpg_utils import to_path

from helpers.pedigree_from_sample_metadata import ext_to_int_sample_map
from reanalysis.utils import get_json_response


MAX_DEPTH: int = 3
PANELAPP_SERVER = 'https://sample-metadata.populationgenomics.org.au'
TEMPLATE = (
    f'{PANELAPP_SERVER}/api/v1/participant/{{dataset}}/individual-metadata-seqr/json'
)


HPO_RE = re.compile(r'HP:[0-9]+')
PANELS_ENDPOINT = 'https://panelapp.agha.umccr.org/api/v1/panels/'


def get_panels(endpoint: str = PANELS_ENDPOINT) -> dict[str, set[int]]:
    """
    query panelapp, and collect panels by HPO term
    returns a dict of {HPO Term: [panel_ID, panel_ID],}
    """
    hpo_dict = defaultdict(set)

    while endpoint:
        endpoint_data = get_json_response(endpoint)
        for panel in endpoint_data['results']:

            # can be split over multiple strings
            relevant_disorders = ' '.join(panel['relevant_disorders'] or [])
            for match in re.findall(HPO_RE, relevant_disorders):
                hpo_dict[match].add(str(panel['id']))

        # cycle through additional pages
        # why don't GEL make the panelapp API public...
        if endpoint_data['next']:
            endpoint = endpoint_data['next']
        else:
            break

    return dict(hpo_dict)


def read_hpo_tree(obo_file: str) -> networkx.MultiDiGraph:
    """
    takes the obo file and creates a graph
    """
    return read_obo(obo_file, ignore_obsolete=False)


def match_hpo_terms(
    panel_map: dict[str, set[int]],
    hpo_tree: networkx.MultiDiGraph,
    hpo_str: str,
    max_layer_delta: int = 3,
    layers_scanned: int = 0,
    selections: set | None = None,
) -> set[int]:
    """
    get a list of panels which are relevant for this HPO
    this includes a manual recursive implementation of the edge traversal
    main reason is to take superseded terms into account

    instead of just checking parent(s), we also check if a term is obsolete
    if so, we instead check each replacement term

    for live terms we recurse on all parents

    this could benefit from some cache-ing, but that will be hard with the
    layer argument

    relevant usage guide:
    https://github.com/dhimmel/obonet/blob/main/examples/go-obonet.ipynb
    """

    if selections is None:
        selections = set()

    # identify identical match and select the panel
    if hpo_str in panel_map:
        selections.update(panel_map[hpo_str])

    if layers_scanned >= max_layer_delta:
        return selections

    # if a node is invalid, recursively call this method for each replacement D:
    # there are simpler ways, just none that are as fun to write
    if not hpo_tree.has_node(hpo_str):
        logging.error(f'HPO term was absent from the tree: {hpo_str}')
        return selections

    hpo_node = hpo_tree.nodes[hpo_str]
    if hpo_node.get('is_obsolete', 'false') == 'true':
        for hpo_term in hpo_node.get('replaced_by', []):
            selections.update(
                match_hpo_terms(
                    panel_map,
                    hpo_tree,
                    hpo_term,
                    max_layer_delta,
                    layers_scanned + 1,
                    selections,
                )
            )
    # search for parent(s), even if the term is obsolete
    for hpo_term in hpo_node.get('is_a', []):
        selections.update(
            match_hpo_terms(
                panel_map,
                hpo_tree,
                hpo_term,
                max_layer_delta,
                layers_scanned + 1,
                selections,
            )
        )
    return selections


def query_and_parse_metadata(dataset_name: str) -> dict:
    """
    takes the seqr metadata and parses out the relevant HPO data
    Parameters
    ----------
    dataset_name : string, the project dataset key to use

    Returns
    -------
    all project metadata, parsed into a dict
    """

    hpo_dict = {}
    seqr_api = SeqrApi()
    participant_meta = seqr_api.get_individual_metadata_for_seqr(project=dataset_name)

    for row in participant_meta['rows']:

        # take the family ID and all HPO terms
        hpo_string = row.get('hpo_terms_present')
        hpo_list = hpo_string.split(',') if hpo_string else []
        hpo_dict[row['individual_id']] = {
            'family_id': row['family_id'],
            'hpo_terms': hpo_list,
        }

    return hpo_dict


def match_hpos_to_panels(
    hpo_to_panel_map: dict,
    hpo_graph: networkx.MultiDiGraph,
    all_hpos: set,
    max_depth: int | None = None,
) -> dict:
    """
    take all the hpo terms, and run each against the panel~hpo matcher

    Parameters
    ----------
    hpo_to_panel_map :
    hpo_graph :
    all_hpos : set of all unique hpo terms
    max_depth : optional overriding graph traversal depth

    Returns
    -------
    a dictionary linking all HPO terms to a corresponding set of Panel IDs
    """

    hpo_to_panels = {}
    for hpo in all_hpos:
        panel_ids = match_hpo_terms(
            panel_map=hpo_to_panel_map,
            hpo_tree=hpo_graph,
            hpo_str=hpo,
            max_layer_delta=max_depth or MAX_DEPTH,
        )
        hpo_to_panels[hpo] = panel_ids

    return hpo_to_panels


def get_unique_hpo_terms(participants_hpo: dict) -> set:
    """
    get all the unique HPO terms across this cohort
    Parameters
    ----------
    participants_hpo :

    Returns
    -------
    all unique HPO terms
    """
    all_hpos = set()
    for participant_dict in participants_hpo.values():
        all_hpos.update(participant_dict['hpo_terms'])
    return all_hpos


def match_participants_to_panels(
    participant_hpos: dict, hpo_panels: dict, participant_map: dict
) -> dict:
    """
    take the two maps of Participants: HPOs, and HPO: Panels
    blend the two to find panels per participant

    For each participant, find any HPO terms which were matched to panels
    for each matched term, add the panel(s) to the participant's private set

    Parameters
    ----------
    participant_hpos :
    hpo_panels :
    participant_map : a lookup of external to CPG ID

    Returns
    -------

    """
    final_dict = {}
    for participant, party_data in participant_hpos.items():
        for participant_key in participant_map.get(participant, [participant]):
            final_dict[participant_key] = {
                'panels': set(),
                'external_id': participant,
                **party_data,
            }
            for hpo_term in party_data['hpo_terms']:
                if hpo_term in hpo_panels:
                    final_dict[participant_key]['panels'].update(hpo_panels[hpo_term])

    return final_dict


def main(dataset: str, output_path: str, obo: str):
    """
    main method linking all component methods together

    Parameters
    ----------
    dataset : the dataset name in metamist
    output_path : path to write output file to
    obo : path to the HPO obo ontology tree
    """

    # get a dictionary of HPO terms to panel IDs
    panels_by_hpo = get_panels()
    hpo_tree = read_hpo_tree(obo_file=obo)

    # pull metadata from metamist/api content
    participants_hpo = query_and_parse_metadata(dataset_name=dataset)

    # obtain a lookup of Ext. ID to CPG ID
    reverse_lookup = ext_to_int_sample_map(project=dataset)

    # mix & match the HPOs, panels, and participants
    # this will be a little complex to remove redundant searches
    # e.g. multiple participants & panels may have the same HPO terms
    # so only search once for each HPO term
    unique_hpos = get_unique_hpo_terms(participants_hpo)
    hpo_to_panels = match_hpos_to_panels(
        hpo_to_panel_map=panels_by_hpo, hpo_graph=hpo_tree, all_hpos=unique_hpos
    )
    participant_panels = match_participants_to_panels(
        participants_hpo, hpo_to_panels, participant_map=reverse_lookup
    )

    with to_path(output_path).open('w') as handle:
        json.dump(participant_panels, handle, indent=4, default=list)


if __name__ == '__main__':

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )
    parser = ArgumentParser()
    parser.add_argument('-d', '--dataset', type=str, help='the dataset to process')
    parser.add_argument('-o', '--output', type=str, help='file path to write output to')
    parser.add_argument('--obo', required=True, help='path to the HPO .obo tree file')
    args = parser.parse_args()
    main(dataset=args.dataset, output_path=args.output, obo=args.obo)
