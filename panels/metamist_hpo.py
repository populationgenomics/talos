"""
script for pulling relevant deets from metamist
"""
import logging
import json
import os
import sys

from argparse import ArgumentParser

import networkx
from sample_metadata.apis import SeqrApi
from sample_metadata.model.export_type import ExportType

from cpg_utils import to_path
from panels.panel_hpo import read_hpo_tree, match_hpo_terms, get_panels

# from reanalysis.utils import get_json_response

MAX_DEPTH: int = 3
PANELAPP_SERVER = 'https://sample-metadata.populationgenomics.org.au'
TEMPLATE = (
    f'{PANELAPP_SERVER}/api/v1/participant/{{dataset}}/individual-metadata-seqr/json'
)


def get_participants(dataset_name: str) -> None:
    """
    This SM-API endpoint isn't currently working, so avoid using this method for now

    Parameters
    ----------
    dataset_name :

    Returns
    -------

    """
    seqr_api = SeqrApi()
    participants = seqr_api.get_individual_metadata_for_seqr(
        project=dataset_name, export_type=ExportType('json')
    )
    return participants


def get_participants_temp(dataset_name: str) -> dict:
    """
    a stopgap method, using cached data
    Parameters
    ----------
    dataset_name :

    Returns
    -------

    """
    metadata_file = os.path.join(os.path.dirname(__file__), f'{dataset_name}-meta.json')
    with open(metadata_file, encoding='utf-8') as handle:
        dictionary = json.load(handle)
    return dictionary


def parse_metadata(participant_meta: dict) -> dict:
    """
    takes the seqr metadata and parses out the relevant HPO data
    Parameters
    ----------
    participant_meta :

    Returns
    -------

    """
    hpo_dict = {}

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
    take all the hpo terms
    Parameters
    ----------
    hpo_to_panel_map :
    hpo_graph :
    all_hpos : set of all unique hpo terms
    max_depth : optional overriding graph traversal depth

    Returns
    -------

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


def match_participants_to_panels(participant_hpos: dict, hpo_panels: dict) -> dict:
    """
    take the two maps of Participants: HPOs, and HPO: Panels
    blend the two to find panels per participant

    For each participant, find any HPO terms which were matched to panels
    for each matched term, add the panel(s) to the participant's private set

    Parameters
    ----------
    participant_hpos :
    hpo_panels :

    Returns
    -------

    """
    final_dict = {}
    for participant, hpo_list in participant_hpos.items():
        final_dict[participant] = set()
        for hpo_term in hpo_list:
            if hpo_term in hpo_panels:
                final_dict[participant].update(hpo_panels[hpo_term])

    return final_dict


def main(dataset: str, output_path: str):
    """

    Parameters
    ----------
    dataset :
    output_path :

    Returns
    -------

    """

    # get a dictionary of HPO terms to panel IDs
    panels_by_hpo = get_panels()
    hpo_tree = read_hpo_tree()

    # pull metadata from metamist/api content
    participant_metadata = get_participants_temp(dataset)
    participants_hpo = parse_metadata(participant_metadata)

    # mix & match the HPOs, panels, and participants
    # this will be a little complex to reduce/remove redundancy
    # e.g. multiple participants & panels may have the same HPO terms
    unique_hpos = get_unique_hpo_terms(participants_hpo)
    hpo_to_panels = match_hpos_to_panels(
        hpo_to_panel_map=panels_by_hpo, hpo_graph=hpo_tree, all_hpos=unique_hpos
    )
    participant_panels = match_participants_to_panels(participants_hpo, hpo_to_panels)

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
    args = parser.parse_args()
    main(dataset=args.dataset, output_path=args.output)
