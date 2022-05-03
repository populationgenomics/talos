"""
generate a ped file on the fly using the sample-metadata api client
optional argument will remove all family associations
    - family structure removal enforces singleton structure during this MVP
optional argument will replace missing parents with '0' to be valid PLINK structure
"""


from typing import Dict, List, Union
from itertools import product
import logging
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


def get_fat_pedigree(
    pedigree_dicts: List[Dict[str, Union[str, List[str]]]],
    sample_to_cpg_dict: Dict[str, List[str]],
    singles: bool,
    plink: bool,
) -> List[Dict[str, List[str]]]:
    """
    swaps sample IDs for the internal CPG values
    optionally make all members unrelated singletons
    optionally substitute missing member IDs for '0'

    :param pedigree_dicts:
    :param sample_to_cpg_dict:
    :param singles:
    :param plink:
    :return:
    """

    new_entries = []

    for counter, ped_entry in enumerate(pedigree_dicts, 1):

        if ped_entry['individual_id'] not in sample_to_cpg_dict:
            continue

        # update the sample IDs
        ped_entry['individual_id'] = sample_to_cpg_dict[ped_entry['individual_id']]

        # remove parents and assign an individual sample ID
        if singles:
            ped_entry['paternal_id'] = ['0' if plink else '']
            ped_entry['maternal_id'] = ['0' if plink else '']
            ped_entry['family_id'] = str(counter)

        else:
            ped_entry['paternal_id'] = sample_to_cpg_dict.get(
                ped_entry['paternal_id'], ['0' if plink else '']
            )
            ped_entry['maternal_id'] = sample_to_cpg_dict.get(
                ped_entry['maternal_id'], ['0' if plink else '']
            )

        new_entries.append(ped_entry)

    return new_entries


def write_fat_pedigree(fat_pedigree: List[Dict[str, List[str]]], output: str):
    """
    take the pedigree data, and write out as a correctly formatted PED file
    this permits the situation where we have multiple possible samples per individual

    :param fat_pedigree: the PED content with sample lists
    :param output: file location to create
    """
    with open(output, 'w', encoding='utf-8') as handle:
        handle.write('\t'.join(PED_KEYS) + '\n')
        for entry in fat_pedigree:
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


def get_pedigree_for_project(project: str) -> List[Dict[str, str]]:
    """
    fetches the project pedigree from sample-metadata
    list, one dict per participant
    :param project:
    :return: all content retrieved from API
    """

    return FamilyApi().get_pedigree(project=project)


def ext_to_int_sample_map(project: str) -> Dict[str, List[str]]:
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

    sample_map = {}
    for (
        participant,
        sample,
    ) in ParticipantApi().get_external_participant_id_to_internal_sample_id(
        project=project
    ):
        sample_map.setdefault(participant, []).append(sample)
    return sample_map


@click.command()
@click.option(
    '--project',
    help='the name of the project to use in API queries',
)
@click.option(
    '--singles',
    default=False,
    is_flag=True,
    help='remake the pedigree as singletons',
)
@click.option(
    '--plink',
    default=False,
    is_flag=True,
    help='make a plink format file',
)
@click.option(
    '--output',
    help='write the new PED file here',
)
def main(project: str, singles: bool, plink: bool, output: str):
    """

    :param project: may be able to retrieve this from the environment
    :param singles: whether to split the pedigree(s) into singletons
    :param plink: whether to write the file as PLINK.fam format
    :param output: path to write new PED file
    """

    # get the list of all pedigree members as list of dictionaries
    logging.info('Pulling all pedigree members')
    pedigree_dicts = get_pedigree_for_project(project=project)

    # endpoint gives list of tuples e.g. [['A1234567_proband', 'CPG12341']]
    # parser returns a dictionary, arbitrary # sample IDs per participant
    logging.info('pulling internal-external sample mapping')
    sample_to_cpg_dict = ext_to_int_sample_map(project=project)

    logging.info('updating pedigree sample IDs to internal')
    fat_pedigree = get_fat_pedigree(
        pedigree_dicts=pedigree_dicts,
        sample_to_cpg_dict=sample_to_cpg_dict,
        singles=singles,
        plink=plink,
    )

    logging.info('writing new PED file to "%s"', output)
    write_fat_pedigree(fat_pedigree, output)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=E1120
