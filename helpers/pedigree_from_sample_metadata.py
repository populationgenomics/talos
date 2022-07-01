"""
generate a ped file on the fly using the sample-metadata api client
optional argument will remove all family associations
    - family structure removal enforces singleton structure during this MVP
optional argument will replace missing parents with '0' to be valid PLINK structure
"""


from collections import defaultdict
from itertools import product
import hashlib
import json
import logging
from typing import Union

import click

from sample_metadata.apis import FamilyApi, ParticipantApi


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
    pedigree_dicts: list[dict[str, Union[str, list[str]]]],
    sample_to_cpg_dict: dict[str, list[str]],
    make_singletons: bool,
    plink_format: bool,
) -> list[dict[str, list[str]]]:
    """
    Take the pedigree entry representations from the pedigree endpoint
    translates sample IDs of all members to CPG values
    sample IDs are replaced with lists, containing all permutations
    where multiple samples exist

    :param pedigree_dicts:
    :param sample_to_cpg_dict:
    :param make_singletons: make all members unrelated singletons
    :param plink_format: substitute missing member IDs for '0'
    :return:
    """

    new_entries = []
    failures: list[str] = []

    for counter, ped_entry in enumerate(pedigree_dicts, 1):

        if ped_entry['individual_id'] not in sample_to_cpg_dict:
            failures.append(ped_entry['individual_id'])

        # update the sample IDs
        ped_entry['individual_id'] = sample_to_cpg_dict[ped_entry['individual_id']]

        # remove parents and assign an individual sample ID
        if make_singletons:
            ped_entry['paternal_id'] = ['0' if plink_format else '']
            ped_entry['maternal_id'] = ['0' if plink_format else '']
            ped_entry['family_id'] = str(counter)

        else:
            ped_entry['paternal_id'] = sample_to_cpg_dict.get(
                ped_entry['paternal_id'], ['0' if plink_format else '']
            )
            ped_entry['maternal_id'] = sample_to_cpg_dict.get(
                ped_entry['maternal_id'], ['0' if plink_format else '']
            )

        new_entries.append(ped_entry)

    if failures:
        raise Exception(
            f'Samples were available from the Pedigree endpoint, '
            f'but no ID translation was available: {",".join(failures)}'
        )
    return new_entries


def write_ped_with_permutations(
    ped_with_permutations: list[dict[str, list[str]]], output: str
):
    """
    take the pedigree data, and write out as a correctly formatted PED file
    this permits the situation where we have multiple possible samples per individual

    :param ped_with_permutations: the PED content with sample lists
    :param output: file location to create
    """
    with open(output, 'w', encoding='utf-8') as handle:
        for entry in ped_with_permutations:
            for sample, mother, father in product(
                entry['individual_id'], entry['paternal_id'], entry['maternal_id']
            ):
                handle.write(
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


def generate_reverse_lookup(mapping_digest: dict[str, list[str]]) -> dict[str, str]:
    """
    :param mapping_digest: created by ext_to_int_sample_map
    :return:
    """

    return {
        sample: participant
        for participant, samples in mapping_digest.items()
        for sample in samples
    }


def hash_reduce_dicts(
    pedigree_dicts: list[dict[str, str]], hash_threshold: int
) -> list[dict[str, str]]:
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


@click.command()
@click.option('--project', help='Project name to use in API queries')
@click.option(
    '--singletons',
    default=False,
    is_flag=True,
    help='remake the pedigree as singletons',
)
@click.option(
    '--plink',
    default=False,
    is_flag=True,
    help='make a plink format file (.fam, .ped is the default)',
)
@click.option('--output', help='prefix for writing all outputs to')
@click.option(
    '--hash_threshold',
    help=(
        'Integer 0-100 representing the % of families to include, e.g. 15'
        'will result in the retention of 15% of families'
    ),
    default=None,
    type=int,
)
def main(
    project: str,
    singletons: bool,
    plink: bool,
    output: str,
    hash_threshold: int | None = None,
):
    """

    :param project: may be able to retrieve this from the environment
    :param singletons: whether to split the pedigree(s) into singletons
    :param plink: whether to write the file as PLINK.fam format
    :param hash_threshold:
    :param output: path to write new PED file
    """

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
        make_singletons=singletons,
        plink_format=plink,
    )

    # store a way of reversing this lookup in future
    reverse_lookup = generate_reverse_lookup(sample_to_cpg_dict)
    with open(f'{output}_external_lookup.json', 'w', encoding='utf-8') as handle:
        json.dump(reverse_lookup, handle, indent=4)

    pedigree_output_path = f'{output}.{"fam" if plink else "ped"}'
    logging.info('writing new PED file to "%s"', pedigree_output_path)
    write_ped_with_permutations(ped_with_permutations, pedigree_output_path)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=E1120
