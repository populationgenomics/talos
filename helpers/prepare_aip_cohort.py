"""
master script for preparing a run

- takes a category; exomes or genomes
- queries for & builds the pedigree
- optionally takes a seqr metadata file
- copies relevant files to GCP
- generates the cohort-specific TOML file
"""


from argparse import ArgumentParser
from collections import defaultdict
from itertools import product
import hashlib
import json
import logging
import os

import toml

from cpg_utils import to_path, Path
from sample_metadata.apis import FamilyApi, ParticipantApi

from reanalysis.utils import read_json_from_path


BUCKET_TEMPLATE = 'gs://cpg-{dataset}-test/reanalysis'
LOCAL_TEMPLATE = 'inputs/{dataset}'


def get_seqr_details(seqr_meta: str, local_root, remote_root) -> tuple[str, str]:
    """
    process the seqr data, write locally, copy remotely
    Args:
        seqr_meta ():
        local_root ():
        remote_root ():

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
    filename = f'seqr_{"exome_" if "exome" in project_id else ""}processed.json'

    local_path = local_root / filename
    remote_path = remote_root / filename
    with local_path.open('w') as handle:
        json.dump(parsed, handle, indent=4, sort_keys=True)
    with remote_path.open('w') as handle:
        json.dump(parsed, handle, indent=4, sort_keys=True)

    return project_id, str(remote_path)


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
    pedigree_dicts: list[dict],
    sample_to_cpg_dict: dict,
) -> list[dict]:
    """
    Take the pedigree entry representations from the pedigree endpoint
    translates sample IDs of all members to CPG values
    sample IDs are replaced with lists, containing all permutations
    where multiple samples exist

    :param pedigree_dicts:
    :param sample_to_cpg_dict:
    :return:
    """

    new_entries = []
    failures: list[str] = []

    for ped_entry in pedigree_dicts:

        if ped_entry['individual_id'] not in sample_to_cpg_dict:
            failures.append(ped_entry['individual_id'])

        # update the sample IDs
        ped_entry['individual_id'] = sample_to_cpg_dict[ped_entry['individual_id']]

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
    ped_with_permutations: list[dict], local_dir: Path, remote_dir: Path
) -> str:
    """
    take the pedigree data, and write out as a correctly formatted PED file
    this permits the situation where we have multiple possible samples per individual

    Args:
        ped_with_permutations (): the PED content with sample lists
        local_dir ():
        remote_dir ():

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

    # write to a local file
    (local_dir / 'pedigree.ped').write_text(single_line)

    # also write to the remote path
    remote_path = remote_dir / 'pedigree.ped'
    remote_path.write_text(single_line)
    logging.info(f'Wrote pedigree with {len(ped_lines)} lines to {remote_path}')

    # pass back the remote file path
    return str(remote_path)


def get_pedigree_for_project(project: str) -> list[dict[str, str]]:
    """
    fetches the project pedigree from sample-metadata
    list, one dict per participant
    :param project:
    :return: all content retrieved from API
    """

    return FamilyApi().get_pedigree(project=project)


def ext_to_int_sample_map(project: str) -> dict[str, list[str]]:
    """
    fetches the participant-sample mapping, so external IDs can be translated
    to the corresponding CPG ID

    This endpoint returns a list of tuples, each containing a participant ID
    and corresponding sample ID

    This originally had an expectation of a 1:1 participant:sample relationship
    that will not hold in general, as numerous samples have multiple samples

    for each participant, create a list of all possible samples

    :param project:
    :return: the mapping dictionary
    """

    sample_map = defaultdict(list)
    for (
        participant,
        sample,
    ) in ParticipantApi().get_external_participant_id_to_internal_sample_id(
        project=project
    ):
        sample_map[participant].append(sample)
    return sample_map


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
    local_lookup = to_path(local_dir) / 'external_lookup.json'
    with local_lookup.open('w') as handle:
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

    :param pedigree_dicts:
    :param hash_threshold: int
    :return:
    """

    reduced_pedigree = []

    for member in pedigree_dicts:
        family_id = member['family_id']
        family_bytes = family_id.encode('utf-8')
        hash_int = int(hashlib.sha1(family_bytes).hexdigest(), 16)
        if hash_int % 100 >= hash_threshold:
            continue
        reduced_pedigree.append(member)

    return reduced_pedigree


def main(project: str, hash_threshold: int = 100, seqr_file: str | None = None):
    """
    Who runs the world? main()
    Returns:

    """

    local_root = to_path(LOCAL_TEMPLATE.format(dataset=project))
    remote_root = to_path(BUCKET_TEMPLATE.format(dataset=project))

    cohort_config = {'dataset_specific': {}}
    if seqr_file:
        project_id, seqr_file = get_seqr_details(seqr_file, local_root, remote_root)
        cohort_config['dataset_specific'] = {
            'seqr_instance': 'https://seqr.populationgenomics.org.au',
            'seqr_project': project_id,
            'seqr_lookup': seqr_file,
        }

    # get the list of all pedigree members as list of dictionaries
    logging.info('Pulling all pedigree members')
    pedigree_dicts = get_pedigree_for_project(project=project)

    # if a threshold is provided, reduce the families present
    if isinstance(hash_threshold, int):
        pedigree_dicts = hash_reduce_dicts(pedigree_dicts, hash_threshold)

    # endpoint gives list of tuples e.g. [['A1234567_proband', 'CPG12341']]
    # parser returns a dictionary, arbitrary # sample IDs per participant
    logging.info('pulling internal-external sample mapping')
    sample_to_cpg_dict = ext_to_int_sample_map(project=project)

    logging.info('updating pedigree sample IDs to internal')
    ped_with_permutations = get_ped_with_permutations(
        pedigree_dicts=pedigree_dicts,
        sample_to_cpg_dict=sample_to_cpg_dict,
    )

    # store a way of reversing this lookup in future
    reverse_lookup = process_reverse_lookup(sample_to_cpg_dict, local_root, remote_root)
    cohort_config['dataset_specific']['external_lookup'] = reverse_lookup

    ped_file = process_pedigree(ped_with_permutations, local_root, remote_root)

    local_config = local_root / 'cohort_config.toml'
    with local_config.open('w') as handle:
        toml.dump(cohort_config, handle)

    print(
        f"""
    --config {local_config}
    --pedigree {ped_file}
    """
    )


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    parser = ArgumentParser()
    parser.add_argument(
        '--project', help='Project name to use in API queries', required=True
    )
    parser.add_argument(
        '--hash_threshold',
        help=(
            'Integer 0-100 representing the percentage of families to '
            'include, e.g. 15 will result in the retention of 15pc of families'
        ),
        type=int,
        default=100,
    )
    parser.add_argument('--seqr', help='optional, seqr JSON file', required=False)
    args = parser.parse_args()
    main(project=args.project, hash_threshold=args.hash_threshold, seqr_file=args.seqr)
