"""
import the HPO obo tree
"""
import os
import re
from collections import defaultdict

import networkx
from obonet import read_obo

from reanalysis.query_panelapp import get_json_response


HPO_RE = re.compile(r'HP:[0-9]+')
OBO_FILE = os.path.join(os.path.dirname(__file__), 'hpo_terms.obo')


def get_panels(endpoint: str) -> dict[str, list[int]]:
    """
    query panelapp, and collect panels by HPO term
    returns a dict of {HPO Term: [panel_ID, panel_ID],}
    """
    hpo_dict = defaultdict(list)

    while endpoint:
        endpoint_data = get_json_response(endpoint)
        for panel in endpoint_data['results']:

            # can be split over multiple strings
            relevant_disorders = ' '.join(panel['relevant_disorders']) or ''
            for match in re.findall(HPO_RE, relevant_disorders):
                hpo_dict[match].append(panel['id'])

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
    return read_obo(obo_file)


def match_hpo_terms(
    panel_map: dict[str, list[int]],
    hpo_tree: networkx.MultiDiGraph,
    hpo_str: str,
    max_layer_delta: int = 3,
) -> list[int]:
    """
    get a list of panels which are relevant for this HPO
    example here is HP:0002194 (a grandchild of HP:0012758, which is mapped)

    relevant usage guide:
    https://github.com/dhimmel/obonet/blob/main/examples/go-obonet.ipynb
    """
    selections = []

    # identical match
    if hpo_str in panel_map:
        print('found exact match')
        selections.extend(panel_map[hpo_str])

    layers_scanned = 0
    # with each iteration, move back a level and check for match on 'parent'
    for _current, parent in networkx.bfs_edges(hpo_tree, hpo_str):
        hpo_str = parent
        layers_scanned += 1
        if layers_scanned == max_layer_delta:
            break
        if hpo_str in panel_map:
            print(f'found match to {hpo_str} at layer {layers_scanned}')
            selections.extend(panel_map[hpo_str])

    return selections


if __name__ == '__main__':
    panels_by_hpo = get_panels('https://panelapp.agha.umccr.org/api/v1/panels/')
    obo_tree = read_hpo_tree(OBO_FILE)

    # get relevant HPO terms from somewhere... metamist?

    # match at 2-layers abstracted
    print(match_hpo_terms(panels_by_hpo, obo_tree, 'HP:0002194'))
    # >>> [250]
    # >>> found match to HP:0012758 at layer 2

    # expecting exact match
    print(match_hpo_terms(panels_by_hpo, obo_tree, 'HP:0012758', max_layer_delta=0))
    # >>> [250]
    # >>> found exact match
