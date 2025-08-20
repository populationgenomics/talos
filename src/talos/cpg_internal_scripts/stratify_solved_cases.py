"""
this will take a file containing family IDs where we know the case
is solved, so it should be removed from future Talos consideration

This is provided as an aggregate across a number of projects, and needs
to be digested down to a per-project object
"""

from sys import argv

import toml
from tenacity import retry

from metamist.graphql import gql, query

PROJECT_QUERY = gql(
    """
    query MyQuery {
        myProjects {
            dataset
        }
    }
    """,
)
ID_QUERY = gql(
    """
    query IdQuery($project: String!) {
      project(name: $project) {
        pedigree
        sequencingGroups {
          id
          sample {
            participant {
              externalId
            }
          }
        }
      }
    }""",
)


@retry()
def get_data_for_project(project: str):
    # find all the samples in each project
    response = query(ID_QUERY, variables={'project': project})

    # find affected individuals
    aff_dict = get_affected_per_family(response['project']['pedigree'])

    # create a mapping of sample ID to participant ID
    sample_map = {sg['sample']['participant']['externalId']: sg['id'] for sg in response['project']['sequencingGroups']}

    solves_project_ids = []

    for solve in solved_fams:
        if solve in aff_dict:
            for ext_id in aff_dict[solve]:
                if ext_id in sample_map:
                    solves_project_ids.append(sample_map[ext_id])
                elif ext_id.replace('_proband', '') in sample_map:
                    solves_project_ids.append(sample_map[ext_id.replace('_proband', '')])

    return sorted(solves_project_ids)


def get_affected_per_family(pedigree: list[dict]):
    """
    given a pedigree, return a dict of family ID to affected individuals

    Args:
        pedigree ():

    Returns:
        a dict of family ID to affected individuals (ext IDs)
    """
    affected_value_in_ped = 2
    affected: dict[str, list[str]] = {}
    for member in pedigree:
        if member['affected'] == affected_value_in_ped:
            affected.setdefault(member['family_id'], []).append(member['individual_id'])
    return affected


# find all my projects
all_projects_of_interest = {dataset['dataset'] for dataset in query(PROJECT_QUERY)['myProjects']}

project_dict = {}
with open(argv[1], encoding='utf-8') as handle:
    solved_fams: set[str] = {x.strip() for x in handle.readlines()}

for project in all_projects_of_interest:
    if 'test' in project or 'training' in project:
        continue

    # write this to a dict
    project_dict[project] = get_data_for_project(project)

with open('output_solves.toml', 'w') as handle:
    toml.dump(project_dict, handle)
