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

from metamist.graphql import gql, query

from cpg_utils import to_path

from reanalysis.utils import get_json_response


PANELAPP_SERVER = 'https://sample-metadata.populationgenomics.org.au'
TEMPLATE = (
    f'{PANELAPP_SERVER}/api/v1/participant/{{dataset}}/individual-metadata-seqr/json'
)

HPO_KEY = 'HPO Terms (present)'
HPO_RE = re.compile(r'HP:[0-9]+')
MAX_DEPTH: int = 3
PANELS_ENDPOINT = 'https://panelapp.agha.umccr.org/api/v1/panels/'


def get_panels(endpoint: str = PANELS_ENDPOINT) -> dict[str, set[int]]:
    """
    query panelapp, and collect panels by HPO term

    Args:
        endpoint (str): URL for panels

    Returns:
        dict: {HPO_Term: [panel_ID, panel_ID],}
    """

    hpo_dict = defaultdict(set)

    while endpoint:
        endpoint_data = get_json_response(endpoint)
        for panel in endpoint_data['results']:
            # can be split over multiple strings
            relevant_disorders = ' '.join(panel['relevant_disorders'] or [])
            for match in re.findall(HPO_RE, relevant_disorders):
                hpo_dict[match].add(int(panel['id']))
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


def query_and_parse_metadata(dataset: str) -> tuple[dict, set[str]]:
    """
    gql query, pull out family details and HPO terms
    may be a little overloaded at the moment
    Args:
        dataset (str): dataset name

    Returns:
        dict of per-participant details, and set of all HPO terms
    """

    query_string = gql(
        """
        query MyQuery($project: String!) {
            project(name: $project) {
                sequencingGroups {
                    sample {
                        participant {
                            phenotypes
                            externalId
                            families {
                                externalId
                            }
                        }
                    }
                    id
                }
            }
        }"""
    )

    result = query(query_string, variables={'project': dataset})
    hpo_dict = {}
    all_hpo: set[str] = set()
    # pylint: disable=unsubscriptable-object
    for sg in result['project']['sequencingGroups']:
        hpos = set(
            HPO_RE.findall(sg['sample']['participant']['phenotypes'].get(HPO_KEY, ''))
        )
        all_hpo.update(hpos)
        hpo_dict[sg['id']] = {
            'hpo_terms': hpos,
            'family_id': sg['sample']['participant']['families'][0]['externalId'],
            'external_id': sg['sample']['participant']['externalId'],
            'panels': {137},  # baseline panel is always mendeliome
        }
    return hpo_dict, all_hpo


def match_hpos_to_panels(
    hpo_to_panel_map: dict,
    obo_file: str,
    all_hpos: set,
    max_depth: int | None = None,
) -> dict:
    """
    take all the hpo terms, and run each against the panel~hpo matcher

    Args:
        hpo_to_panel_map ():
        obo_file (str): file containing HPO tree
        all_hpos (set[str]): set of all unique hpo terms
        max_depth (int): optional overriding graph traversal depth

    Returns:
        a dictionary linking all HPO terms to a corresponding set of Panel IDs
    """

    hpo_graph = read_obo(obo_file, ignore_obsolete=False)

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


def match_participants_to_panels(participant_hpos: dict, hpo_panels: dict):
    """
    take the two maps of Participants: HPOs, and HPO: Panels
    blend the two to find panels per participant

    For each participant, find any HPO terms which were matched to panels
    for each matched term, add the panel(s) to the participant's private set

    Args:
        participant_hpos ():
        hpo_panels ():
        participant_map ():
    """
    for party_data in participant_hpos.values():
        for hpo_term in party_data['hpo_terms']:
            if hpo_term in hpo_panels:
                party_data['panels'].update(hpo_panels[hpo_term])


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

    # pull metadata from metamist/api content
    participants_hpo, unique_hpos = query_and_parse_metadata(dataset=dataset)

    # mix & match the HPOs, panels, and participants
    # this will be a little complex to remove redundant searches
    # e.g. multiple participants & panels may have the same HPO terms
    # so only search once for each HPO term
    hpo_to_panels = match_hpos_to_panels(
        hpo_to_panel_map=panels_by_hpo, obo_file=obo, all_hpos=unique_hpos
    )
    match_participants_to_panels(participants_hpo, hpo_to_panels)

    with to_path(output_path).open('w') as handle:
        json.dump(participants_hpo, handle, indent=4, default=list)


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
