"""
import the HPO obo tree
parse the panelapp
"""

import logging
import os
import re
import sys

from collections import defaultdict

import networkx
from obonet import read_obo

from reanalysis.query_panelapp import get_json_response


HPO_RE = re.compile(r'HP:[0-9]+')
OBO_FILE = os.path.join(os.path.dirname(__file__), 'hpo_terms.obo')
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
                hpo_dict[match].add(panel['id'])

        # cycle through additional pages
        # why don't GEL make the panelapp API public...
        if endpoint_data['next']:
            endpoint = endpoint_data['next']
        else:
            break

    return dict(hpo_dict)


def read_hpo_tree(obo_file: str = OBO_FILE) -> networkx.MultiDiGraph:
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
        return selections
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


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )
    panels_by_hpo = get_panels()
    obo_tree = read_hpo_tree()

    # get relevant HPO terms from somewhere... metamist?

    # match at 2-layers abstracted
    logging.info(match_hpo_terms(panels_by_hpo, obo_tree, 'HP:0002194'))
    # >>> [250]
    # >>> found match to HP:0012758 at layer 2

    # expecting exact match
    logging.info(
        match_hpo_terms(panels_by_hpo, obo_tree, 'HP:0012758', max_layer_delta=0)
    )
    # >>> [250]
    # >>> found exact match
