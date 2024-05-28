#! /usr/bin/env python3

"""
This script generates a loaded PED file for use in AIP
This headerless output file is a TSV in the format:
    Family ID
    Individual ID
    Paternal ID (or 0 is missing)
    Maternal ID (or 0 is missing)
    Sex (1=male; 2=female; other=unknown)
    Phenotype (-9/0=missing, 1=unaffected, 2=affected)
    Ext ID (or a repeat of Individual ID)
    HPO terms, one per column, no limit on number of columns

This is the CPG-specific implementation of this process
"""

import re
from argparse import ArgumentParser

from metamist.graphql import gql, query

HPO_KEY = 'HPO Terms (present)'
HPO_RE = re.compile(r'HP:[0-9]+')
PARTICIPANT_QUERY = gql(
    """
query MyQuery($project: String!) {
  project(name: $project) {
    pedigree
    sequencingGroups {
      id
      sample {
        participant {
          externalId
          phenotypes
        }
      }
    }
  }
}""",
)

# result = query(PARTICIPANT_QUERY, variables={'project': dataset})


def get_data_from_metamist(project: str) -> list[list[str]]:
    """
    Query metamist for the required data
    Args:
        project ():

    Returns:
        returns the new Ped contents, ready to be written to a file
        each row is a list of Strings, forming a regular pedigree
        followed by an external ID (can be a repeat of internal ID)
        then an arbitrary number of columns for HPO terms
    """

    # to store the new PED entries
    ped_entries: list[list[str]] = []

    # first get a lookup of Int IDs to Ext IDs
    result = query(PARTICIPANT_QUERY, variables={'project': project})

    # will map CPG IDs to external project participant IDs
    int_to_ext: dict[str, str] = {}

    # will map each CPG ID to its set of HPO terms (can be empty)
    cpg_to_hpos: dict[str, set[str]] = {}

    # iterate over all SG Entities
    for sg in result['project']['sequencingGroups']:

        int_to_ext[sg['id']] = sg['sample']['participant']['externalId']
        cpg_to_hpos[sg['id']] = set(HPO_RE.findall(sg['sample']['participant']['phenotypes'].get(HPO_KEY, '')))

    # iterate over the pedigree entities, forming a list from each. List elements:
    # Family ID, Individual ID, Paternal, Maternal, Sex, Affection status, Ext ID (or repeat), HPOs (or not if absent)
    for entry in result['project']['pedigree']:
        # eg = {
        #     "family_id": "19W000139",
        #     "individual_id": "19W000139",
        #     "paternal_id": None,
        #     "maternal_id": None,
        #     "sex": 1,
        #     "affected": 2,
        #     "notes": "",
        # }
        ped_row: list[str] = [
            entry['family_id'],
            entry['individual_id'],
            entry['paternal_id'] or '0',
            entry['maternal_id'] or '0',
            str(entry['sex']),
            str(entry['affected']),
            int_to_ext.get(entry['individual_id'], entry['individual_id']),  # defaults to internal ID again
        ]

        # if there are recorded HPOs, extend the row with them
        if hpos := cpg_to_hpos.get(entry['individual_id']):
            ped_row.extend(sorted(hpos))

        ped_entries.append(ped_row)

    return ped_entries


if __name__ == '__main__':
    parser = ArgumentParser(description='Generate a PED file for AIP')
    parser.add_argument('dataset', help='The dataset to query for')
    parser.add_argument('output', help='The output file')
    args = parser.parse_args()

    new_ped_rows = get_data_from_metamist(args.dataset)

    # write a headless TSV file, as an extended PED format
    with open(args.output, 'w', encoding='utf-8') as handle:
        for line in new_ped_rows:
            handle.write('\t'.join(line) + '\n')
