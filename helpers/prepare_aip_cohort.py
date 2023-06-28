"""
master script for preparing a run

- takes a category; exomes or genomes
- queries for & builds the pedigree
- optionally takes a seqr metadata file
- copies relevant files to GCP
- generates the cohort-specific TOML file
- tweaks for making singleton versions of the given cohort
"""

from argparse import ArgumentParser
from itertools import product
from typing import Any
import hashlib
import json
import logging
import os
import re
import toml

from cpg_utils import to_path, Path
from cpg_utils.config import get_config

from metamist.graphql import gql, query

from reanalysis.utils import read_json_from_path
from helpers.hpo_panel_matching import (
    get_panels,
    get_unique_hpo_terms,
    match_hpos_to_panels,
    query_and_parse_metadata,
)
from helpers.utils import ext_to_int_sample_map

BUCKET_TEMPLATE = 'gs://cpg-{dataset}-test-analysis/reanalysis'
LOCAL_TEMPLATE = 'inputs/{dataset}'
OBO_DEFAULT = os.path.join(os.path.dirname(__file__), 'hpo_terms.obo')

MAX_DEPTH: int = 3

HPO_RE = re.compile(r'HP:[0-9]+')
PANELS_ENDPOINT = 'https://panelapp.agha.umccr.org/api/v1/panels/'
PRE_PANEL_PATH = to_path(__file__).parent.parent / 'reanalysis' / 'pre_panelapp.json'


def match_participants_to_panels(
    participant_hpos: dict, hpo_panels: dict, participant_map: dict
) -> dict:
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

    final_dict = {}
    for participant, party_data in participant_hpos.items():
        for participant_key in participant_map.get(participant, [participant]):
            final_dict[participant_key] = {
                'panels': {137},  # always default to mendeliome
                'external_id': participant,
                **party_data,
            }
            for hpo_term in party_data['hpo_terms']:
                if hpo_term in hpo_panels:
                    final_dict[participant_key]['panels'].update(hpo_panels[hpo_term])

    # now populate the missing samples?
    for ext_id, int_ids in participant_map.items():
        if ext_id not in participant_hpos:
            for each_id in int_ids:
                final_dict[each_id] = {
                    'panels': {137},  # always default to mendeliome
                    'external_id': ext_id,
                    'hpo_terms': [],
                }

    return final_dict


def get_seqr_details(
    seqr_meta: str, local_root, remote_root, exome: bool = False
) -> tuple[str, str]:
    """
    process the seqr data, write locally, copy remotely
    Args:
        seqr_meta ():
        local_root ():
        remote_root ():
        exome ():

    Returns:

    """

    if not os.path.exists(seqr_meta):
        raise FileExistsError(f'Input file {seqr_meta} inaccessible')

    details_dict = read_json_from_path(seqr_meta)

    # map CPG ID to individual GUID
    parsed = {
        sample['sampleId']: sample['familyGuid']
        for seqr_sample_id, sample in details_dict['samplesByGuid'].items()
    }

    project_id = {
        sample['projectGuid'] for sample in details_dict['samplesByGuid'].values()
    }
    assert len(project_id) == 1, f'Multiple projects identified: {project_id}'
    project_id = project_id.pop()

    logging.info(f'{len(parsed)} families in seqr metadata')
    filename = f'seqr_{"exome_" if exome else ""}processed.json'

    with (local_root / filename).open('w') as handle:
        json.dump(parsed, handle, indent=4, sort_keys=True)
    with (remote_root / filename).open('w') as handle:
        json.dump(parsed, handle, indent=4, sort_keys=True)

    return project_id, str(remote_root / filename)


# the keys provided by the SM API, in the order to write in output
PED_KEYS = [
    'family_id',
    'individual_id',
    'paternal_id',
    'maternal_id',
    'sex',
    'affected',
]


def get_ped_with_permutations(
    pedigree_dicts: list[dict], sample_to_cpg_dict: dict, make_singletons: bool = False
) -> list[dict]:
    """
    Take the pedigree entry representations from the pedigree endpoint
    translates sample IDs of all members to CPG values
    sample IDs are replaced with lists, containing all permutations
    where multiple samples exist
    Optionally, overwrite all family structures and render as singletons
    If we're running singletons, remove all unaffected

    Args:
        pedigree_dicts ():
        sample_to_cpg_dict ():
        make_singletons ():

    Returns:
        list of dicts representing rows of the pedigree
    """

    new_entries = []
    failures: list[str] = []

    # enumerate to get ints - use these as family IDs if singletons
    for counter, ped_entry in enumerate(pedigree_dicts, 1):
        if ped_entry['individual_id'] not in sample_to_cpg_dict:
            failures.append(ped_entry['individual_id'])

        # update the sample IDs
        ped_entry['individual_id'] = sample_to_cpg_dict[ped_entry['individual_id']]

        if make_singletons:
            # skip unaffected singletons
            if ped_entry['affected'] == 0:
                continue
            ped_entry['paternal_id'] = ['0']
            ped_entry['maternal_id'] = ['0']
            ped_entry['family_id'] = str(counter)
        else:
            # remove parents and assign an individual sample ID
            ped_entry['paternal_id'] = sample_to_cpg_dict.get(
                ped_entry['paternal_id'], ['0']
            )
            ped_entry['maternal_id'] = sample_to_cpg_dict.get(
                ped_entry['maternal_id'], ['0']
            )

        new_entries.append(ped_entry)

    if failures:
        logging.error(
            f'Samples were available from the Pedigree endpoint, '
            f'but no ID translation was available: {",".join(failures)}'
        )
    return new_entries


def process_pedigree(
    ped_with_permutations: list[dict],
    local_dir: Path,
    remote_dir: Path,
    singletons: bool = False,
) -> str:
    """
    take the pedigree data, and write out as a correctly formatted PED file
    this permits the situation where we have multiple possible samples per individual

    Args:
        ped_with_permutations (): the PED content with sample lists
        local_dir ():
        remote_dir ():
        singletons (bool):

    Returns:
        The path we're writing the remote data to
    """

    ped_lines = []
    for entry in ped_with_permutations:
        for sample, mother, father in product(
            entry['individual_id'], entry['paternal_id'], entry['maternal_id']
        ):
            ped_lines.append(
                '\t'.join(
                    [
                        entry['family_id'],
                        sample,
                        mother,
                        father,
                        str(entry['sex']),
                        str(entry['affected']),
                    ]
                )
                + '\n'
            )

    # condense all lines into one
    single_line = ''.join(ped_lines)

    # make the pedigree name
    ped_name = f'{"singletons_" if singletons else ""}pedigree.ped'

    # write to a local file
    (local_dir / ped_name).write_text(single_line)

    # also write to the remote path
    remote_path = remote_dir / ped_name
    remote_path.write_text(single_line)
    logging.info(f'Wrote pedigree with {len(ped_lines)} lines to {remote_path}')

    # pass back the remote file path
    return str(remote_path)


def get_pedigree_for_project(project: str) -> list[dict[str, str]]:
    """
    fetches the project pedigree from sample-metadata
    list, one dict per participant

    Args:
        project (str): project/dataset to use in query

    Returns:
        All API returned content
    """
    ped_query = gql(
        """
    query MyQuery($project: String!) {
        project(name: $project) {
            pedigree
        }
    }
    """
    )
    # pylint: disable=unsubscriptable-object
    response: dict[str, Any] = query(ped_query, variables={'project': project})
    return response['project']['pedigree']


def process_reverse_lookup(
    mapping_digest: dict[str, list], local_dir: Path, remote_dir: Path
) -> str:
    """

    Args:
        mapping_digest ():
        local_dir ():
        remote_dir ():

    Returns:

    """

    clean_dict = {
        sample: participant
        for participant, samples in mapping_digest.items()
        for sample in samples
    }
    with (to_path(local_dir) / 'external_lookup.json').open('w') as handle:
        json.dump(clean_dict, handle, indent=4)

    remote_lookup = to_path(remote_dir) / 'external_lookup.json'
    with remote_lookup.open('w') as handle:
        json.dump(clean_dict, handle, indent=4)

    return str(remote_lookup)


def hash_reduce_dicts(pedigree_dicts: list[dict], hash_threshold: int) -> list[dict]:
    """
    hashes the family ID of each member of the Pedigree
    Normalises the Hash value to the range 0 - 99
    if the normalised value exceeds the threshold, remove

    Args:
        pedigree_dicts ():
        hash_threshold ():

    Returns:

    """

    logging.info(f'Reducing families to {hash_threshold}%')

    reduced_pedigree = []

    for member in pedigree_dicts:
        family_id = member['family_id']
        family_bytes = family_id.encode('utf-8')
        hash_int = int(hashlib.sha1(family_bytes).hexdigest(), 16)
        if hash_int % 100 >= hash_threshold:
            continue
        reduced_pedigree.append(member)

    return reduced_pedigree


def main(
    project: str,
    obo: str,
    seqr_file: str | None = None,
    exome: bool = False,
    singletons: bool = False,
):
    """
    Who runs the world? main()

    Args:
        project (str): project to query for
        obo (str): path to an HPO graph file
        seqr_file (str):
        exome ():
        singletons ():
    """

    local_root = to_path(LOCAL_TEMPLATE.format(dataset=project))

    if not local_root.exists():
        local_root.mkdir(parents=True, exist_ok=True)

    remote_root = to_path(BUCKET_TEMPLATE.format(dataset=project))

    # get the list of all pedigree members as list of dictionaries
    logging.info('Pulling all pedigree members')
    pedigree_dicts = get_pedigree_for_project(project=project)

    # if a threshold is provided, reduce the families present
    hash_threshold = get_config()['cohorts'][project].get('cohort_percentage', 100)
    if hash_threshold != 100:
        pedigree_dicts = hash_reduce_dicts(pedigree_dicts, hash_threshold)

    # endpoint gives list of tuples e.g. [['A1234567_proband', 'CPG12341']]
    # parser returns a dictionary, arbitrary # sample IDs per participant
    logging.info('pulling internal-external sample mapping')
    sample_to_cpg_dict = ext_to_int_sample_map(project=project)

    logging.info('updating pedigree sample IDs to internal')
    ped_with_permutations = get_ped_with_permutations(
        pedigree_dicts=pedigree_dicts,
        sample_to_cpg_dict=sample_to_cpg_dict,
        make_singletons=singletons,
    )

    # store a way of reversing this lookup in future
    reverse_lookup = process_reverse_lookup(sample_to_cpg_dict, local_root, remote_root)

    # try a more universal way of preparing output paths
    path_prefixes = []
    if exome:
        path_prefixes.append('exomes')
    if singletons:
        path_prefixes.append('singleton')

    logging.info(f'Output Prefix:\n---\nreanalysis/{"/".join(path_prefixes)}\n---')

    cohort_config = {
        'dataset_specific': {
            'historic_results': str(
                remote_root / '/'.join(path_prefixes + ['historic_results'])
            ),
            'external_lookup': reverse_lookup,
        }
    }
    if seqr_file:
        project_id, seqr_file = get_seqr_details(
            seqr_file, local_root, remote_root, exome
        )
        cohort_config['dataset_specific'].update(
            {
                'seqr_instance': 'https://seqr.populationgenomics.org.au',
                'seqr_project': project_id,
                'seqr_lookup': seqr_file,
            }
        )

    ped_file = process_pedigree(
        ped_with_permutations, local_root, remote_root, singletons=singletons
    )

    path_prefixes.append('cohort_config.toml')
    cohort_path = local_root / '_'.join(path_prefixes)

    with cohort_path.open('w') as handle:
        toml.dump(cohort_config, handle)
        logging.info(f'Wrote cohort config to {cohort_path}')

    # pull metadata from metamist/api content
    participants_hpo = query_and_parse_metadata(dataset_name=project)

    # mix & match the HPOs, panels, and participants
    # this will be a little complex to remove redundant searches
    # e.g. multiple participants & panels may have the same HPO terms
    # so only search once for each HPO term
    hpo_to_panels = match_hpos_to_panels(
        hpo_to_panel_map=get_panels(),
        obo_file=obo,
        all_hpos=get_unique_hpo_terms(participants_hpo),
    )
    participant_panels = match_participants_to_panels(
        participants_hpo, hpo_to_panels, participant_map=sample_to_cpg_dict
    )
    panel_local = local_root / 'participant_panels.json'
    with panel_local.open('w') as handle:
        json.dump(participant_panels, handle, indent=4, default=list)
        logging.info(f'Wrote panel file to {panel_local}')

    panel_remote = remote_root / 'participant_panels.json'
    with panel_remote.open('w') as handle:
        json.dump(participant_panels, handle, indent=4, default=list)
        logging.info(f'Wrote panel file to {panel_remote}')

    # finally, copy the pre-panelapp content if it didn't already exist
    if 'pre_panelapp' in (
        prior := get_config()['cohorts'][project].get('gene_prior', 'MISSING')
    ):
        pre_panelapp = read_json_from_path(PRE_PANEL_PATH)
        with to_path(prior).open('w') as handle:
            json.dump(pre_panelapp, handle, indent=4)
            logging.info(f'Wrote VCGS gene prior file to {prior}')

    logging.info(f'--pedigree {ped_file}')


if __name__ == '__main__':
    logging.basicConfig(level=logging.WARN)
    parser = ArgumentParser()
    parser.add_argument(
        '--project', help='Project name to use in API queries', required=True
    )
    parser.add_argument('--seqr', help='optional, seqr JSON file', required=False)
    parser.add_argument(
        '--obo', default=OBO_DEFAULT, help='path to the HPO .obo tree file'
    )
    parser.add_argument('-e', help='cohort is exomes', action='store_true')
    parser.add_argument(
        '--singletons', help='cohort as singletons', action='store_true'
    )
    args = parser.parse_args()
    main(
        project=args.project,
        obo=args.obo,
        seqr_file=args.seqr,
        exome=args.e,
        singletons=args.singletons,
    )
