"""
utils to prevent circular imports
"""

from collections import defaultdict
from typing import Any

from metamist.graphql import gql, query


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

    mapping_query = gql(
        """
    query MyQuery($project: String!) {
        project(name: $project) {
            participants {
                samples {
                    sequencingGroups {
                        id
                    }
                }
                externalId
            }
        }
    }
    """
    )
    response: dict[str, Any] = query(mapping_query, variables={'project': project})

    sample_map = defaultdict(list)
    for participant_block in response['project']['participants']:
        ext_id = participant_block['externalId']
        samples = [
            sam['id']
            for sg in participant_block['samples']
            for sam in sg['sequencingGroups']
        ]
        sample_map[ext_id].extend(samples)
    return sample_map
