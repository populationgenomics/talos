#!/usr/bin/env python3

"""
Revised process for generating pedigree and phenotype data for Talos.

- Go to Metamist
- Query for all participants in a project, and their phenotypes
- Query for the pedigree data, and replace with family external IDs
- Integrate the data into a 7-column pedigree format, with ','-delimited HPO terms in the last column
"""

import re
from argparse import ArgumentParser

from cloudpathlib.anypath import to_anypath
from loguru import logger

from metamist.graphql import gql, query

HPO_KEY = 'HPO Terms (present)'
HPO_RE = re.compile(r'HP:\d{7}')

PARTICIPANT_QUERY = gql(
    """
query Participants($project: String!, $sequencing_type: String!, $technology: String!) {
  project(name: $project) {
    sequencingGroups(
        technology: {eq: $technology},
        type:  {eq: $sequencing_type},
        activeOnly: {eq: true}
    ) {
      id
      sample {
        participant {
          externalIds
          phenotypes
        }
      }
    }
  }
}""",
)

# whole-dataset pedigree, with no restriction on sequencing type or technology
PEDIGREE_QUERY = gql(
    """
query Pedigree($project: String!) {
  project(
      name: $project
  ) {
    pedigree(
      replaceWithFamilyExternalIds: true
    )
  }
}""",
)


def get_pedigree_data(dataset: str) -> dict[str, dict]:
    """
    Query for the Metamist pedigree data for a given dataset. Inclusive of all sequencing types and technologies.
    Example data format:
    {
      "family_id": "family_id",
      "individual_id": "participant_id",
      "paternal_id": null,
      "maternal_id": null,
      "sex": 2,
      "affected": 1,
      "notes": ""
    },
    """
    result = query(PEDIGREE_QUERY, variables={'project': dataset})
    return {entry['individual_id']: entry for entry in result['project']['pedigree']}


def get_participant_data(
    dataset: str,
    seq_type: str,
    tech: str = 'short-read',
) -> dict[str, dict[str, str | list[str]]]:
    """
    Query for the Metamist participant data for a given dataset.

    Output format:
    {
        'participant_id': {
            'sgid': 'sequencing_group_id',
            'phenotypes': ['HP:0000001', 'HP:0000002', ...],
        },
        ...
    }
    """
    result = query(
        PARTICIPANT_QUERY,
        variables={'project': dataset, 'sequencing_type': seq_type, 'technology': tech},
    )

    # index this stuff back on Participant ID, so we can match it to the pedigree data
    participant_data = {}
    for sg in result['project']['sequencingGroups']:
        # default external ID is the empty-string key on the externalIds dict... IKR, I didn't implement this.
        participant_id = sg['sample']['participant']['externalIds']['']
        if not participant_id:
            logger.warning(f'Sequence group {sg["id"]} has no participant ID, skipping')
            continue

        # grab the HPO terms, if they exist. Defaults to an empty list if not absent/empty
        hpos = set(HPO_RE.findall(sg['sample']['participant']['phenotypes'].get(HPO_KEY, '')))

        participant_data[participant_id] = {
            'sgid': sg['id'],
            'phenotypes': hpos,
        }

    return participant_data


def integrate_participant_data(
    pedigree_data: dict[str, dict],
    participant_data: dict[str, dict[str, str | list[str]]],
) -> list[str]:
    """Integrate the participant data into the pedigree data, building up a 7-column pedigree format."""

    # index all SG IDs by participant ID for quick access
    sg_id_map = {participant_id: details['sgid'] for participant_id, details in participant_data.items()}

    pedigree_lines = []

    # iterate over the participant data, which is the more restrictive query (by sequencing type and technology)
    for participant_id, participant_details in participant_data.items():
        if participant_id not in pedigree_data:
            logger.warning(f'Participant ID {participant_id} not found in pedigree data, skipping')
            continue

        # get this participant's pedigree data
        participant_pedigree = pedigree_data[participant_id]

        sg_id = participant_details['sgid']

        # get the participant's phenotypes as a comman-separated string
        phenotypes = participant_details['phenotypes']
        hpo_terms = ','.join(sorted(phenotypes)) if phenotypes else ''

        # create a new entry with the HPO terms and mother/father IDs (subbed out for 0s if not present)
        new_entry = [
            pedigree_data[participant_id]['family_id'],
            sg_id,
            sg_id_map.get(participant_pedigree['maternal_id'], '0'),
            sg_id_map.get(participant_pedigree['paternal_id'], '0'),
            str(participant_pedigree['sex']),
            str(participant_pedigree['affected']),
            hpo_terms,
        ]

        pedigree_lines.append('\t'.join(new_entry) + '\n')

    return pedigree_lines


def main(output: str, dataset: str, seq_type: str, tech: str = 'short-read'):
    """Assemble a HPO-annotated Pedigree from Metamist."""

    pedigree_data = get_pedigree_data(dataset=dataset)
    participant_data = get_participant_data(dataset=dataset, seq_type=seq_type, tech=tech)
    pedigree_lines = integrate_participant_data(
        pedigree_data=pedigree_data,
        participant_data=participant_data,
    )

    with to_anypath(output).open('w', encoding='utf-8') as handle:
        for line in pedigree_lines:
            handle.write(line)


def cli_main():
    parser = ArgumentParser(description='Generate a PED file for Talos')
    parser.add_argument('--dataset', help='The dataset to query for')
    parser.add_argument('--output', help='The output file')
    parser.add_argument('--type', help='Sequencing type (exome or genome)', required=True)
    parser.add_argument('--tech', help='Sequencing technology', default='short-read')
    args = parser.parse_args()

    main(output=args.output, dataset=args.dataset, seq_type=args.type, tech=args.tech)


if __name__ == '__main__':
    cli_main()
