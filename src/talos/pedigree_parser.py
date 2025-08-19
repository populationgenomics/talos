"""
A general purpose Pedigree parser, developed for Talos to facilitate validation of an extended-format Pedigree file.
This is modelled on a classic 6-column Pedigree file format, with the following specification:

- Columns are Tab-separated
- A couple of standard 'null' values are used: "0", "-", or empty string
- Column 1 is the Family ID, used to connect members of a family, and present family ID in the results, mandatory
- Column 2 is the Sample ID, must match exactly with the sample ID in the VCFs/Variant Data, mandatory
- Column 3 is the Father ID, which must match the Sample ID of a participant in the same family, or be empty
- Column 4 is the Mother ID, which must match the Sample ID of a participant in the same family, or be empty
- Column 5 is the Sex, represented as "1" for Male, "2" for Female, or "0" for Unknown
- Column 6 is the Affected Status, represented as a 1 or null value for Unaffected, or "2" for Affected
- Column 7 is optional, and contains HPO terms for the affected participant
    - "," or ";"-delimited HPO terms, which will be parsed into a set of HPO terms for the affected participant, e.g.
        - "HP:0000118,HP:0001250" for a participant with two HPO terms
        - "HP:0001250" for a participant with one HPO term
"""

import re
from dataclasses import dataclass

from cloudpathlib.anypath import to_anypath
from loguru import logger

HPO_RE = re.compile(r'HP:\d{7}')

NULL_INT = 0
NULL_PARTICIPANT = '0'
MISSING_ID = {'0', 0}
VALID_SEX = {'1', '2', '0'}
MALE_SEX = {'male', '1'}
FEMALE_SEX = {'female', '2'}
UNKNOWN_SEX = {'0', 'unknown', 'other'}
NULL_VALUES = {'0', 0, 'other', '-', '', 'unaffected'}
VALID_AFFECTED = {'0', '1', '2'}
GRUDGINGLY_VALID_AFFECTED = {'affected', 'true'}

EXPECTED_NUM_COLUMNS = 7
AFFECTED_NUM = FEMALE_SEX_INT = 2
UNAFFECTED_NUM = MALE_SEX_INT = 1


@dataclass
class Participant:
    """Contain all the information about a participant in the pedigree."""

    family_id: str
    sample_id: str
    father_id: str | None
    mother_id: str | None
    sex: int
    affected: int
    hpo_terms: set[str]

    def __str__(self):
        """String representation of the participant, used when writing a Pedigree."""
        return f'{self.family_id}\t{self.sample_id}\t{self.father_id}\t{self.mother_id}\t{self.sex}\t{self.affected}'

    @property
    def is_affected(self):
        """Boolean property to check if the participant is affected."""
        return self.affected == AFFECTED_NUM

    @property
    def is_not_affected(self):
        """Boolean property to check if the participant is not affected."""
        return self.affected != AFFECTED_NUM

    @property
    def is_male(self):
        """Boolean property to check if the participant is male."""
        return self.sex == MALE_SEX_INT

    @property
    def is_female(self):
        """Boolean property to check if the participant is female."""
        return self.sex == FEMALE_SEX_INT


PEDIGREE_DATA = dict[str, Participant]


class PedigreeParser:
    def __init__(self, pedigree_path: str):
        """
        A Pedigree object, because the world needs another Pedigree-parsing class.
        This takes a single argument, a path to a pedigree file, and parses it into a dictionary of Participant objects.

        This module attempts to do a very lenient validation of the pedigree file, allowing some common substitutions.
        If there are errors in the pedigree file, they will be logged, and a ValueError will be raised.

        If no errors are found, a PedigreeParser instance will have two main attributes:
          - self.participants: a dictionary of sample IDs: Participants, representing all participants in the pedigree
          - self.by_family: a dictionary of family IDs: list of Participants, grouping all participants by family

        class methods can be used to retrieve select data:
          - get_affected_members: returns an {ID: Participant} dictionary of only affected participants
          - get_affected_member_ids: returns a set of sample IDs of affected participants
          - as_singletons: returns a dictionary of participants as singletons, stripping parental information
          - strip_pedigree_to_samples: reduces self.participants to only the requested sample IDs
          - write_pedigree: writes the participants to a PLINK format 6-column PED file

        Args:
            pedigree_path (str): Path to the pedigree file, which is expected to be a 6/7-column tab-separated file
        """
        self.pedigree_path = pedigree_path
        self.parsing_issues: list[str] = []
        self.participants: PEDIGREE_DATA = self.read_pedigree()
        self.by_family: dict[str, list[Participant]] = self.get_participants_by_family()

        if self.parsing_issues:
            logger.warning('Issues found during pedigree parsing:')
            for issue in self.parsing_issues:
                logger.warning(issue)
            raise ValueError('Errors found during pedigree parsing, see log for details')

        if len(self.participants) == 0:
            raise ValueError('No valid participants found in the pedigree file!')

    def get_participants_by_family(self) -> dict[str, list[Participant]]:
        """Index the pre-collected participants by family ID."""
        by_family: dict[str, list[Participant]] = {}
        for participant in self.participants.values():
            by_family.setdefault(participant.family_id, []).append(participant)
        return by_family

    def get_affected_members(self) -> PEDIGREE_DATA:
        """
        Return only the affected participants from the provided dictionary of participants.

        This will return all participants with an affected status of 1
        """
        return {
            sample_id: participant for sample_id, participant in self.participants.items() if participant.is_affected
        }

    def get_affected_member_ids(self) -> set[str]:
        """Return a set of sample IDs for all affected participants in the pedigree."""
        return set(self.get_affected_members().keys())

    def get_all_sample_ids(self) -> set[str]:
        """Return a set of all sample IDs in the pedigree."""
        return set(self.participants.keys())

    def as_singletons(self) -> PEDIGREE_DATA:
        """
        Return a dictionary of participants as singletons, i.e. each participant is a separate entry in the dictionary.

        This will return all participants, stripped of parental information
        """
        return_participants = {}
        for sample_id, participant in self.participants.items():
            # create a new instance from the original details, don't risk overwriting the original Participant
            return_participants[sample_id] = Participant(
                family_id=participant.family_id,
                sample_id=participant.sample_id,
                father_id=NULL_PARTICIPANT,
                mother_id=NULL_PARTICIPANT,
                sex=participant.sex,
                affected=participant.affected,
                hpo_terms=participant.hpo_terms,
            )
        return return_participants

    def set_participants(self, participants: PEDIGREE_DATA) -> None:
        """
        Set the participants attribute to a new dictionary of participants.
        This is used to replace the participants with a new dictionary, e.g. after subsetting or creating singletons.
        """
        self.participants = participants
        self.by_family = self.get_participants_by_family()

    def strip_pedigree_to_samples(self, only_participants: list[str] | set[str]) -> PEDIGREE_DATA:
        """
        Takes the Pedigree data, and strips it back to only samples in the `only_participants` object
        Returns this as a new object. If required this can overwrite the main data using self.set_participants().
        A situation for this is when a provided pedigree is a superset of the participants we have variant data for. By
        stripping the pedigree down to only the samples we've really seen, we can use all pedigree participants with no
        need to constantly re-check.
        """
        return {key: value for key, value in self.participants.items() if key in only_participants}

    def write_pedigree(
        self,
        output_path: str,
        participants: PEDIGREE_DATA | None = None,
        only_participants: list[str] | None = None,
    ) -> None:
        """
        Write the participants to a PLINK-format Fam/PED file at `output_path`.

        The output will be a 6-column PED file, with the following columns:
        - Family ID
        - Sample ID
        - Father ID
        - Mother ID
        - Sex
        - Affected Status

        the only_participants parameter can be used to limit the output to a specific set of participants.
        """

        logger.info(f'Writing PED file to: {output_path}')

        # if a participant dictionary is not provided, use the class's participants
        if participants is None:
            participants = self.participants

        with to_anypath(output_path).open('w') as filehandle:
            for sample_id, participant in participants.items():
                if only_participants and sample_id not in only_participants:
                    continue  # skip participants not in the only_participants list

                # write the participant to the file using the __str__ method
                filehandle.write(str(participant) + '\n')

        logger.info(f'PED file written to: {output_path}')

    def read_pedigree(self, prune_missing_parents: bool = True) -> PEDIGREE_DATA:
        """
        Read a PED file and return a dictionary of Participant objects.

        The default validation action here is to strip off any mother/father IDs if the mother/father are not present in
        the pedigree file. If prune_missing_parents is set to False, this will instead cause a failure.
        """

        participants = {}

        logger.info(f'Processing PED file: {self.pedigree_path}')

        with to_anypath(self.pedigree_path).open() as filehandle:
            for line in filehandle:
                if line.startswith('#') or not line.strip():
                    continue  # Skip comments and empty lines

                parts = line.rstrip().split('\t')

                if len(parts) > EXPECTED_NUM_COLUMNS:
                    logger.error(f'Skipping malformed line: {line.strip()}')
                    continue

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
                if len(parts) == EXPECTED_NUM_COLUMNS:
                    phenotypes_block = parts[6].rstrip()
                    # we ask for , as a delimiter, but we detect semicolon as well
                    delimiter = ';' if ';' in phenotypes_block else ','
                    candidate_hpo_terms = {term.strip() for term in phenotypes_block.split(delimiter)}

                    for each_hpo in candidate_hpo_terms:
                        if HPO_RE.match(each_hpo):
                            hpo_terms.add(each_hpo)
                        else:
                            logger.warning(f'Invalid HPO term found: {each_hpo} in line: {line.strip()}')

                sex_int = self.validate_sex(sex_str=sex, sample_id=sample_id)
                affected_int = self.validate_affected(aff_str=affected, sample_id=sample_id)

                # Create a Participant object, even if it contained some validation issues
                participants[sample_id] = Participant(
                    family_id=family_id,
                    sample_id=sample_id,
                    father_id=father_id if father_id not in MISSING_ID else '0',
                    mother_id=mother_id if mother_id not in MISSING_ID else '0',
                    sex=sex_int,
                    affected=affected_int,
                    hpo_terms=hpo_terms,
                )

        # Now validate parents, and validate that all listed IDs are present in the participants
        # we're making this valus a String consistently, but we allow '0' as a valid value for missing parents
        for each_id, each_participant in participants.items():
            if each_participant.father_id not in participants and each_participant.father_id != '0':
                message = f'Participant {each_id} has an invalid father ID: {each_participant.father_id}'
                if prune_missing_parents:
                    logger.info(message)
                    each_participant.father_id = '0'
                else:
                    self.parsing_issues.append(message)

            # same check on the mother
            if each_participant.mother_id not in participants and each_participant.mother_id != '0':
                message = f'Participant {each_id} has an invalid mother ID: {each_participant.mother_id}'
                if prune_missing_parents:
                    logger.info(message)
                    each_participant.mother_id = '0'
                else:
                    self.parsing_issues.append(message)

        return participants

    def validate_sex(self, sex_str: str, sample_id: str) -> int:
        """Validate the provided value for Sex."""

        if sex_str in VALID_SEX:
            return int(sex_str)

        lower_sex = sex_str.lower()
        if lower_sex in MALE_SEX:
            return MALE_SEX_INT
        if lower_sex in FEMALE_SEX:
            return FEMALE_SEX_INT
        if lower_sex in UNKNOWN_SEX:
            return NULL_INT

        # a value is returned here to prevent the model validation from failing, but this will be caught and reported
        self.parsing_issues.append(f'Invalid Sex provided! Sample {sample_id}: {sex_str}')
        return NULL_INT

    def validate_affected(self, aff_str: str, sample_id: str) -> int:
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
            return UNAFFECTED_NUM

        if lower_aff in GRUDGINGLY_VALID_AFFECTED:
            # this is a warning, but we still return 1 to indicate affected status
            logger.warning(
                f'Grudgingly valid Affected status provided, please correct data! Sample {sample_id}: {aff_str}',
            )
            return AFFECTED_NUM

        self.parsing_issues.append(f'Invalid Affected status provided! Sample {sample_id}: {aff_str}')
        return NULL_INT
