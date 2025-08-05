"""
Pedigree parsing in Talos has a limited scope, so a simple parser has been implemented here

- Each affected participant will be considered separately, i.e. two affected siblings are considered as two trios
- For each Proband, we want to connect to parents to apply MOI testing to the immediate nuclear family

This is modelled on a classic 6-column Pedigree file format, with the following specification:

- Columns are Tab-separated
- A couple of standard 'null' values are used: "0", "-", or empty string
- Column 1 is the Family ID, used to connect members of a family, and present family ID in the results, mandatory
- Column 2 is the Sample ID, must match exactly with the sample ID in the VCFs/Variant Data, mandatory
- Column 3 is the Father ID, which must match the Sample ID of a participant in the same family, or be empty
- Column 4 is the Mother ID, which must match the Sample ID of a participant in the same family, or be empty
- Column 5 is the Sex, represented as "1" for Male, "2" for Female, or "0" for Unknown
- Column 6 is the Affected Status, represented as a 1 or null value for Unaffected, or one of the values below:
    - "2" for Affected, matching a classic pedigree format
    - ";"-delimited HPO terms, which will be parsed into a set of HPO terms for the affected participant, e.g.
        - "HP:0000118;HP:0001250" for a participant with two HPO terms
        - "HP:0001250" for a participant with one HPO term
"""

import re
import sys
from dataclasses import dataclass

from cloudpathlib.anypath import AnyPath, to_anypath
from loguru import logger

from talos import models


HPO_RE = re.compile(r'HP:\d{7}')

NULL_PARTICIPANT = '0'
MISSING_ID = {'0', 0}
VALID_SEX = {'1', '2', '0'}
MALE_SEX = {'male', '1'}
FEMALE_SEX = {'female', '2'}
UNKNOWN_SEX = {'0', 'unknown', 'other'}
NULL_VALUES = {'0', 0, 'other', '-', '', 'unaffected'}
VALID_AFFECTED = {'0', '1', '2'}
GRUDGINGLY_VALID_AFFECTED = {'affected', 'true'}

# contain all issues found during parsing, print once, and comprehensively
ISSUES: list[str] = []

# todo exclude_absent_members(Pedigree, set[str]) -> Pedigree:


# move this to pydantic?
@dataclass
class Participant:
    """Contain all the information about a participant in the pedigree."""

    family_id: str
    sample_id: str
    father_id: str | None
    mother_id: str | None
    sex: int
    affected: int
    hpo_terms: set[str] | None = None

    def __str__(self):
        """String representation of the participant, used when writing a Pedigree."""
        return f'{self.family_id}\t{self.sample_id}\t{self.father_id}\t{self.mother_id}\t{self.sex}\t{self.affected}'


def validate_sex(sex_str: str, sample_id: str) -> int:
    """Validate the provided value for Sex."""

    if sex_str in VALID_SEX:
        return int(sex_str)

    lower_sex = sex_str.lower()
    if lower_sex in MALE_SEX:
        return 1
    if lower_sex in FEMALE_SEX:
        return 2
    if lower_sex in UNKNOWN_SEX:
        return 0

    global ISSUES
    ISSUES.append(f'Invalid Sex provided! Sample {sample_id}: {sex_str}')

    # a value is returned here to prevent the model validation from failing, but this will be caught and reported
    return 0


def validate_affected(aff_str: str, sample_id: str) -> int:
    """
    Validate the provided value for Affected status.
    Return values are 0, 1, or 2, where:
    - 0 means invalid Affected status or unknown
    - 1 means Unaffected
    - 2 means Affected
    """

    if aff_str in VALID_AFFECTED:
        return int(aff_str)

    lower_aff = aff_str.lower()
    # happy to parse these
    if lower_aff in NULL_VALUES:
        return 1

    if lower_aff in GRUDGINGLY_VALID_AFFECTED:
        # this is a warning, but we still return 1 to indicate affected status
        logger.warning(f'Grudgingly valid Affected status provided, please correct data! Sample {sample_id}: {aff_str}')
        return 2

    global ISSUES
    ISSUES.append(f'Invalid Affected status provided! Sample {sample_id}: {aff_str}')
    return 0


def affected_members(participants: dict[str, Participant]) -> dict[str, Participant]:
    """
    Return a only the affected participants from the provided dictionary of participants.

    This will return all participants with an affected status of 1
    """
    return {sample_id: participant for sample_id, participant in participants.items() if participant.affected == 2}


def as_singletons(participants: dict[str, Participant]) -> dict[str, Participant]:
    """
    Return a dictionary of participants as singletons, i.e. each participant is a separate entry in the dictionary.

    This will return all affected participants, stripped of parental information
    """
    return_participants = {}
    for sample_id, participant in participants.items():
        participant.father_id = NULL_PARTICIPANT
        participant.mother_id = NULL_PARTICIPANT
        return_participants[sample_id] = participant
    return return_participants


def write_pedigree(participants: dict[str, Participant], output_path: str, only_participants: list[str] = None) -> None:
    """
    Write the participants to a PED file.

    The output will be a 6-column PED file, with the following columns:
    - Family ID
    - Sample ID
    - Father ID
    - Mother ID
    - Sex
    - Affected Status

    the only_participants parameter can be used to limit the output to a specific set of participants.
    """
    output_path = to_anypath(output_path)

    logger.info(f'Writing PED file to: {output_path}')

    with output_path.open('w') as filehandle:
        for sample_id, participant in participants.items():
            if only_participants and sample_id not in only_participants:
                continue  # skip participants not in the only_participants list

            # write the participant to the file using the __str__ method
            filehandle.write(str(participant) + '\n')

    logger.info(f'PED file written to: {output_path}')


def read_pedigree(filepath: str, prune_missing_parents: bool = True) -> dict[str, Participant]:
    """
    Read a PED file and return a dictionary of Participant objects.

    The default validation action here is to strip off any mother/father IDs if the mother/father are not present in the
    pedigree file. If prune_missing_parents is set to False, this will instead cause a failure.
    """

    participants = {}

    logger.info(f'Processing PED file: {filepath}')

    with to_anypath(filepath).open() as filehandle:
        for line in filehandle:
            if line.startswith('#') or not line.strip():
                continue  # Skip comments and empty lines

            parts = line.rstrip().split('\t')

            if len(parts) > 7:
                logger.error(f'Skipping malformed line: {line.strip()}')
                continue  # Skip malformed lines

            hpo_terms: set[str] = set()

            # extract the mandatory components of the line
            (
                family_id,
                sample_id,
                father_id,
                mother_id,
                sex,
                affected,
            ) = parts[:6]

            # validate the mandatory components
            if any(MISSING_ID.__contains__(each_id) for each_id in [family_id, sample_id]):
                logger.error(f'Skipping participant with missing family or sample ID: {line.strip()}')
                continue

            # we'll wrap around to validate parents once we've gathered all participants

            # gather HPO terms
            if len(parts) == 7:
                phenotypes_block = parts[6].strip()
                # we ask for ; as a delimited, but we detect commas
                delimiter = ';' if ';' in phenotypes_block else ','
                candidate_hpo_terms = {term.strip() for term in phenotypes_block.split(delimiter)}

                for each_hpo in candidate_hpo_terms:
                    if HPO_RE.match(each_hpo):
                        hpo_terms.add(each_hpo)
                    else:
                        logger.warning(f'Invalid HPO term found: {each_hpo} in line: {line.strip()}')

            sex_int = validate_sex(sex_str=sex, sample_id=sample_id)
            affected_int = validate_sex(sex_str=sex, sample_id=sample_id)

            # Create a Participant object, even if it contained some validation issues
            participants[sample_id] = models.PedigreeMember(
                family=family_id,
                id=sample_id,
                father=father_id if father_id not in MISSING_ID else '0',
                mother=mother_id if mother_id not in MISSING_ID else '0',
                sex=sex_int,
                affected=affected_int,
                hpo_terms=hpo_terms,
            )

    # Now validate parents, and validate that all listed IDs are present in the participants
    global ISSUES
    for each_id, each_participant in participants.items():
        if each_participant.father not in participants and each_participant.father != '0':
            message = f'Participant {each_id} has an invalid father ID: {each_participant.father}'
            if prune_missing_parents:
                logger.info(message)
                each_participant.father = 0
            else:
                ISSUES.append(message)

        # same check on the mother
        if each_participant.mother not in participants and each_participant.mother != '0':
            message = f'Participant {each_id} has an invalid mother ID: {each_participant.mother}'
            if prune_missing_parents:
                logger.info(message)
                each_participant.mother = 0
            else:
                ISSUES.append(message)

    if ISSUES:
        logger.warning('Issues found during pedigree parsing:')
        for issue in ISSUES:
            logger.warning(issue)
            sys.exit(1)

    return participants
