"""
master script for preparing a run, can be run from local

- takes a category; exomes or genomes
- queries for & builds the pedigree
- optionally takes a seqr metadata file
- copies relevant files to GCP
- generates the cohort-specific TOML file
- tweaks for making singleton versions of the given cohort
"""
# mypy: ignore-errors
from argparse import ArgumentParser
from itertools import product
import json
import logging
import os
import toml

from cpg_utils import to_path, Path

from metamist.graphql import gql, query

from reanalysis.utils import read_json_from_path, get_cohort_config
from reanalysis.hpo_panel_match import main as hpo_match

BUCKET_TEMPLATE = 'gs://cpg-{dataset}-test-analysis/reanalysis'
LOCAL_TEMPLATE = 'inputs/{dataset}'
OBO_DEFAULT = os.path.join(os.path.dirname(__file__), 'hpo_terms.obo')
PRE_PANEL_PATH = to_path(__file__).parent.parent / 'reanalysis' / 'pre_panelapp.json'

PED_QUERY = gql(
    """
    query PedAndSGs($project: String!) {
        project(name: $project) {
            pedigree
            sequencingGroups(activeOnly: {eq: true}) {
                id
                sample {
                    participant {
                        externalId
                    }
                }
            }
        }
    }
    """
)


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
    one_project_id: str = project_id.pop()

    logging.info(f'{len(parsed)} families in seqr metadata')
    filename = f'seqr_{"exome_" if exome else ""}processed.json'

    with (local_root / filename).open('w') as handle:
        json.dump(parsed, handle, indent=4, sort_keys=True)
    with (remote_root / filename).open('w') as handle:
        json.dump(parsed, handle, indent=4, sort_keys=True)

    return one_project_id, str(remote_root / filename)


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
    pedigree_dicts: list[dict], ext_lookup: dict, make_singletons: bool = False
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
        ext_lookup (dict[str, str]):
        make_singletons ():

    Returns:
        list of dicts representing rows of the pedigree
    """

    new_entries = []
    failures: list[str] = []

    # enumerate to get ints - use these as family IDs if singletons
    for counter, ped_entry in enumerate(pedigree_dicts, 1):
        if ped_entry['individual_id'] not in ext_lookup:
            failures.append(ped_entry['individual_id'])
            continue

        # update the sample IDs
        ped_entry['individual_id'] = ext_lookup[ped_entry['individual_id']]

        if make_singletons:
            # skip unaffected singletons
            if ped_entry['affected'] == 0:
                continue
            ped_entry['paternal_id'] = ['0']
            ped_entry['maternal_id'] = ['0']
            ped_entry['family_id'] = str(counter)
        else:
            # remove parents and assign an individual sample ID
            ped_entry['paternal_id'] = ext_lookup.get(ped_entry['paternal_id'], ['0'])
            ped_entry['maternal_id'] = ext_lookup.get(ped_entry['maternal_id'], ['0'])

        new_entries.append(ped_entry)

    if failures:
        logging.error(
            f'Samples were available from the Pedigree endpoint, '
            f'but no ID translation was available: {", ".join(failures)}'
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
            entry['individual_id'],
            entry['paternal_id'] or ['0'],
            entry['maternal_id'] or ['0'],
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


def get_pedigree_for_project(
    project: str,
) -> tuple[list[dict[str, str]], dict[str, str]]:
    """
    fetches the project pedigree from sample-metadata
    list, one dict per participant

    Args:
        project (str): project/dataset to use in query

    Returns:
        All API returned content
    """
    response = query(PED_QUERY, variables={'project': project})
    pedigree = response['project']['pedigree']
    lookup = {
        sg['sample']['participant']['externalId']: [sg['id']]
        for sg in response['project']['sequencingGroups']
    }
    return pedigree, lookup


def main(
    project: str,
    obo: str,
    seqr_file: str | None = None,
    exome_or_genome: str = 'genome',
    singletons: bool = False,
):
    """
    Who runs the world? main()

    Args:
        project (str): project to query for
        obo (str): path to an HPO graph file
        seqr_file (str):
        exome_or_genome (str): flag to switch exome/genome behaviour
        singletons ():
    """

    local_root = to_path(LOCAL_TEMPLATE.format(dataset=project))

    if not local_root.exists():
        local_root.mkdir(parents=True, exist_ok=True)

    remote_root = to_path(BUCKET_TEMPLATE.format(dataset=project))

    # find HPO-matched panels for each participant
    hpo_panel_dict = hpo_match(
        dataset=project,
        hpo_file=obo,
        panel_out=str(local_root / 'participant_panels.json'),
    )

    panel_remote = remote_root / 'participant_panels.json'
    with panel_remote.open('w') as handle:
        json.dump(hpo_panel_dict, handle, indent=4, default=list)
        logging.info(f'Wrote panel file to {panel_remote}')

    # get the list of all pedigree members as list of dictionaries
    logging.info('Pulling all pedigree members')
    pedigree_dicts, ext_lookup = get_pedigree_for_project(project=project)

    # endpoint gives list of tuples e.g. [['A1234567_proband', 'CPG12341']]
    # parser returns a dictionary, arbitrary # sample IDs per participant
    logging.info('pulling internal-external sample mapping')

    logging.info('updating pedigree sample IDs to internal')
    ped_with_permutations = get_ped_with_permutations(
        pedigree_dicts=pedigree_dicts,
        ext_lookup=ext_lookup,
        make_singletons=singletons,
    )

    # store a way of reversing this lookup in future
    path_prefixes = []
    if exome_or_genome == 'exome':
        path_prefixes.append('exomes')
    if singletons:
        path_prefixes.append('singleton')

    logging.info(f'Output Prefix:\n---\nreanalysis/{"/".join(path_prefixes)}\n---')

    if seqr_file:
        project_id, seqr_file = get_seqr_details(
            seqr_file, local_root, remote_root, exome_or_genome == 'exome'
        )
        seqr_files = {
            'seqr_instance': 'https://seqr.populationgenomics.org.au',
            'seqr_project': project_id,
            'seqr_lookup': seqr_file,
        }
    else:
        seqr_files = {}

    cohort_config = {
        'cohorts': {
            project: {
                exome_or_genome: {
                    'historic_results': str(
                        remote_root / '/'.join(path_prefixes + ['historic_results'])
                    )
                }
                | seqr_files
            }
        },
        'workflow': {'sequencing_type': exome_or_genome},
    }

    ped_file = process_pedigree(
        ped_with_permutations, local_root, remote_root, singletons=singletons
    )

    path_prefixes.append('cohort_config.toml')

    with (local_root / '_'.join(path_prefixes)).open('w') as handle:
        toml.dump(cohort_config, handle)
        logging.info(f'Wrote cohort config to {local_root / "_".join(path_prefixes)}')

    # finally, copy the pre-panelapp content if it didn't already exist
    if 'pre_panelapp' in (
        prior := get_cohort_config(project).get('gene_prior', 'MISSING')
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
    parser.add_argument('--seqr', help='optional, seqr JSON file')
    parser.add_argument(
        '--obo', default=OBO_DEFAULT, help='path to the HPO .obo tree file'
    )
    parser.add_argument('-e', help='cohort is exomes', action='store_true')
    parser.add_argument(
        '--singletons', help='cohort as singletons', action='store_true'
    )
    args = parser.parse_args()
    E_OR_G = 'exome' if args.e else 'genome'
    main(
        project=args.project,
        obo=args.obo,
        seqr_file=args.seqr,
        exome_or_genome=E_OR_G,
        singletons=args.singletons,
    )
