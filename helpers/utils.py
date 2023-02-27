"""
utils to prevent circular imports
"""


from collections import defaultdict

from sample_metadata.apis import ParticipantApi


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
