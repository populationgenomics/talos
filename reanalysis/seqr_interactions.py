"""
api experiments with an auth-blocked seqr instance

this is just a sketch - no fancy requests handling
"""

import json
from collections import defaultdict

import click
import requests

import google.auth.exceptions
import google.auth.transport.requests
from google.oauth2 import service_account

from cpg_utils.config import get_config
from sample_metadata.apis import ProjectApi

from reanalysis.html_builder import make_coord_string


# grab this from config in future
# SEQR_INSTANCE = 'https://seqr-reanalysis-dev.populationgenomics.org.au'
SEQR_INSTANCE = 'http://localhost:3000'

SEQR_AUDIENCE = (
    '1021400127367-4vch8s8kc9opeg4v14b2n70se55jpao4.apps.googleusercontent.com'
)
CREATE_ENDPOINT = '/api/saved_variant/create'
SAMPLES_ENDPOINT = '/api/project/{projectGuid}/samples/sa'
VARIANTS_ENDPOINT = '/api/project/{projectGuid}/saved_variants/'
AIP_VARIANTS_ENDPOINT = '/api/project/{projectGuid}/aip_saved_variants'

# mapping of AIP categories to seqr tag names
AIP_CATEGORY_FLAGS = {
    'category_1': 'AIP Cat 1',
    'category_2': 'AIP Cat 2',
    'category_3': 'AIP Cat 3',
    'category_5': 'AIP Cat 5',
    'category_support': 'AIP Cat Support',
}
AIP_DE_NOVO_FLAG = {'category_4': 'AIP Cat 4'}


def get_token(
    cred_file: str = '/Users/mwelland/Downloads/seqr-308602-713cfeb8eeb6.json',
):
    """
    borrowed/stolen from github.com/populationgenomics/sample-metadata/blob/
    81e7dddead12729e78280eaf0911a2cd59e34c1c/scripts/sync_seqr.py#L549
    """

    with open(cred_file, encoding='utf-8') as f:

        info = json.load(f)
        credentials_content = (info.get('type') == 'service_account') and info or None
        credentials = service_account.IDTokenCredentials.from_service_account_info(
            credentials_content, target_audience=SEQR_AUDIENCE
        )
        auth_req = google.auth.transport.requests.Request()
        credentials.refresh(auth_req)
        return credentials.token


def get_project_from_metamist() -> str | None:
    """
    get the seqr project ID for this project/seq type
    Returns:
        str: the Seqr project GUID or None
    """
    project_dataset = get_config()['workflow']['dataset']
    proj_seqr_key = f'seqr-project-{get_config()["workflow"]["sequencing_type"]}'

    # filter the returned projects
    projects = [
        proj
        for proj in ProjectApi().get_seqr_projects()
        if proj.get('name') == project_dataset
    ]

    # WALRUS
    if results := len(projects) != 1:
        raise Exception(f'Expected one project for {project_dataset}, found {results}')

    project = projects[0]

    if proj_seqr_key in project['meta'] and project['meta'].get('is_seqr'):
        return project['meta'].get(proj_seqr_key)

    return None


def get_current_variants_canonical(project: str, header: dict[str, str]) -> dict[dict]:
    """
    queries seqr for the current variants in this project

    Returns:
        dictionary digest -
        FamilyID: {
            VARIANT_STRING: [
                list of tags applied
            ]
        }
    """

    # create a container for these variant details
    tagged_variants = defaultdict(dict)

    seqr_variant_url = f'{SEQR_INSTANCE}{VARIANTS_ENDPOINT.format(projectGuid=project)}'

    saved_var_request = requests.get(seqr_variant_url, headers=header, timeout=60)
    saved_var_request.raise_for_status()

    # collect all the different tag names
    # happy to assume this key exists, otherwise all this is pointless
    tag_lookup = {
        key: tag['name']
        for key, tag in saved_var_request.json()['variantTagsByGuid'].items()
    }

    # the find all prev. tagged variants
    # keep the guid in case we need to update instead of add
    for guid, variant in saved_var_request.json()['savedVariantsByGuid'].items():
        var_id = variant['variantId']
        tags = [tag_lookup[tag] for tag in variant['tagGuids']]

        for family in variant['familyGuids']:
            tagged_variants[family][var_id] = {
                'tags': tags,
                'guid': guid,
            }

    return dict(tagged_variants)


def get_current_variants_aip_redux(project: str, header: dict[str, str]) -> dict[dict]:
    """
    queries seqr for the current variants in this project

    Returns:
        dictionary digest -
        FamilyID: {
            VARIANT_STRING: [
                list of tags applied
            ]
        }
    """

    # create a container for these variant details
    tagged_variants = defaultdict(dict)

    # use the private endpoint
    seqr_variant_url = (
        f'{SEQR_INSTANCE}{AIP_VARIANTS_ENDPOINT.format(projectGuid=project)}'
    )

    saved_var_request = requests.get(seqr_variant_url, headers=header, timeout=60)
    saved_var_request.raise_for_status()

    # the find all prev. tagged variants
    # keep the guid in case we need to update instead of add
    for guid, variant in saved_var_request.json().items():
        var_id = variant['variantId']

        # possibility of no tag names...
        # probably just from me fiddling around
        if not variant['tagNames']:
            continue

        for family in variant['families']:
            tagged_variants[family][var_id] = {
                'tags': set(variant['tagNames']),
                'guid': guid,
            }

    return dict(tagged_variants)


def get_family_ids_for_cpg(project: str, seqr_header: dict) -> dict[str, str]:
    """
    hits the sa endpoint for sample ID mapping
    Args:
        project (): seqr project ID
        seqr_header (): HEADERS

    Returns:
        mapping of cpg to seqr families
    """

    sa_url = f'{SEQR_INSTANCE}{SAMPLES_ENDPOINT.format(projectGuid=project)}'
    get_resp = requests.get(sa_url, headers=seqr_header, timeout=10)
    get_resp.raise_for_status()
    resp_samples = get_resp.json()['samplesByGuid']

    # reverse this later?
    return {sample['sampleId']: sample['familyGuid'] for sample in resp_samples}


def harvest_aip_results(aip_json_path: str, sample_map: dict) -> dict:
    """
    Pull out relevant details from the AIP JSON
    CPG IDs are at the individual level
    Seqr IDs are at the family level
    Here we need to be aware of that relationship, and
        allow for the same variant multiple times against
        the seqr family, and force a set across collected
        flags

    Args:
        aip_json_path ():
        sample_map ():

    Returns:
        SeqrID: {
            variantID: {
                Tag name,
                Tag name,
            },
        },
    """

    aip_holder = defaultdict(dict)

    with open(aip_json_path, encoding='utf-8') as handle:
        aip_results = json.load(handle)

    for cpg_sample, variants in aip_results.items():
        seqr_family = sample_map.get(cpg_sample)
        if seqr_family is None:
            print(f'{cpg_sample} not present in this seqr project')
            continue

        for variant in variants:
            var_id = make_coord_string(variant)

            # core boolean categories
            var_cats = {
                name for cat, name in AIP_CATEGORY_FLAGS.items() if variant[cat]
            }

            # only valid if this sample had a de novo
            for cat, name in AIP_DE_NOVO_FLAG.items():
                if cpg_sample in variant[cat]:
                    var_cats.add(name)

            if var_id in aip_holder[seqr_family]:
                aip_holder[seqr_family][var_id].update(var_cats)

            else:
                aip_holder[seqr_family][var_id] = var_cats

    return aip_holder


def filter_aip_to_new_flags(aip_digest: dict, extant_flags: dict) -> list:
    """
    use the aip results and existing flags to find what new
    ones we need to add
    Args:
        aip_digest ():
        extant_flags ():

    Returns:
        a list of the new flags to add
    """

    flags_to_add = []

    for family_id, variants in aip_digest.items():
        # add all entries as new flags
        if family_id not in extant_flags:
            for var_id, tags in variants.items():
                for tag_name in tags:
                    flags_to_add.append(
                        {'family': family_id, 'tag': tag_name, 'variant': var_id}
                    )
            continue

        # now check to remove redundancy
        existing = extant_flags[family_id]

        for variant_id, tags in variants.items():
            for tag in tags - existing.get(variant_id, set()):
                flags_to_add.append(
                    {'family': family_id, 'tag': tag, 'variant': variant_id}
                )
            continue

    return flags_to_add


def post_flags(flag_list: list[dict], header: dict):
    """
    ping off those new flags
    Args:
        flag_list ():
        header ():

    Returns:

    """

    create_flag_url = f'{SEQR_INSTANCE}{CREATE_ENDPOINT}'
    for new_flag in flag_list:

        var_id = new_flag['variant']
        chrom, pos, ref, alt = var_id.split('-')
        chrom = chrom if 'chr' in chrom else f'chr{chrom}'

        # sending a variant body with 'note' will prevent the tag being generated
        # sending a tag with an existing 'tagGuid' will update an existing object
        # sending a variant body with tags: [{'name': 'BLAH'}]
        # where BLAH is a name (not a guid) of a tag, will generate that tag
        variant_body = {
            'familyGuid': new_flag['family'],
            'tags': [{'name': new_flag['tag'], 'metadata': 'Added by AIP'}],
            'variant': [
                {
                    'alt': alt,
                    'chrom': chrom,
                    'genomeVersion': '38',
                    'pos': int(pos),
                    'ref': ref,
                    'variantId': new_flag['variant'],
                }
            ],
        }

        # create URL isn't parametrized by project
        req = requests.post(
            create_flag_url, headers=header, json=variant_body, timeout=10
        )
        req.raise_for_status()
        if req.ok:
            print(f'Added new tag for {new_flag}')


@click.command()
@click.option('--aip_json')
def main(aip_json: str):
    """
    run all the things!
    Args:
        aip_json (): path to the results of an AIP run

    Returns:

    """

    # get the seqr project GUID for this dataset
    # R0001_validation on my local instance, not in metamist
    my_proj = get_project_from_metamist()
    if my_proj is None:
        raise Exception('No project identified')

    # location of the auth key could be the deployment/container default
    # or passed manually as an argument
    seqr_header = {'Authorization': f'Bearer {get_token()}'}

    # mapping of sample IDs
    sample_map = get_family_ids_for_cpg(project=my_proj, seqr_header=seqr_header)

    # pull all current variants from seqr
    existing_flags = get_current_variants_aip_redux(project=my_proj, header=seqr_header)

    # load AIP results, and back filter to remove already tagged variants
    aip_digest = harvest_aip_results(aip_json_path=aip_json, sample_map=sample_map)

    # back filter to find only the remaining variants
    flags_to_post = filter_aip_to_new_flags(
        aip_digest=aip_digest, extant_flags=existing_flags
    )

    post_flags(flags_to_post, seqr_header)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
