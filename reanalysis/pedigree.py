"""
class implementing bidirectional pedigree methods

1. Get the raw participant data, indexed by sample ID
2. Cast the pedigree members as objects,
    including immediate relationships upwards
    i.e. if the participant has parents, make those objects and
    assign them as the mother/father of this node
3. Assign children to members where appropriate
"""

from typing import Dict, List, Optional, Set, Type
from dataclasses import dataclass
from csv import DictReader
from cloudpathlib import AnyPath


@dataclass
class PedEntry:
    """
    class representing each participant from the PED
    """

    def __init__(self, fam, sample, father, mother, female, affected):
        """

        :param fam:
        :param sample:
        :param father:
        :param mother:
        :param female:
        :param affected:
        """
        self.family: str = fam
        self.sample_id: str = sample

        # parental IDs as Strings or None
        self.father: str = father if father != '' else None
        self.mother: str = mother if mother != '' else None

        self.is_female = female == '2'
        self.affected = affected == '2'


@dataclass
class Participant:
    """
    dataclass representing a person within a family
    Type['Participant'] has to be used in order to have
    a self-referential class...
    """

    details: PedEntry
    mother: Optional[Type['Participant']]
    father: Optional[Type['Participant']]
    children: List[Type['Participant']]
    affected_parents: Set[str]
    unaffected_parents: Set[str]


PED_KEYS = [
    '#Family ID',
    'Individual ID',
    'Paternal ID',
    'Maternal ID',
    'Sex',
    'Affected',
]


class PedigreeParser:
    """
    takes a PED file, and reads into a collection of linked-list like objects
    """

    def __init__(self, pedfile: str):
        """

        :param pedfile: path to a PED file
        """
        self.ped_dict = self.read_participants(pedfile)
        self.participants: Dict[str, Participant] = {}
        for sample_id in self.ped_dict.keys():
            self.populate_participants(sample_id=sample_id)
        self.apply_children()

        # no need for this object, but small memory footprint so who cares
        # del self.ped_dict

    @staticmethod
    def read_participants(ped_file: str) -> Dict[str, PedEntry]:
        """
        Iterates through a PED file, and parses into a dict of Participants
        the dict is indexed by the Participant Sample ID string

        :param ped_file:
        :return:
        """

        with open(AnyPath(ped_file), 'r', encoding='utf-8') as handle:
            ped_reader = DictReader(handle, delimiter='\t')
            participants = {
                party_line[PED_KEYS[1]]: PedEntry(
                    *[party_line.get(key) for key in PED_KEYS]
                )
                for party_line in ped_reader
            }

        return participants

    def populate_participants(self, sample_id: str):
        """
        take a sample ID, and works backwards through the pedigree
        if we find a parent not already made into an object, create
        finally create a Participant for _this_ sample, including
        references to their parents, also as Participant objects
        :param sample_id:
        :return:
        """

        if sample_id in self.participants:
            return

        # create two sets for this participant - all parent sample_ids which are
        # affected and unaffected. This prevents recalculating this list during every
        # MOI test later on
        affected_parents = set()
        unaffected_parents = set()

        ped_sample = self.ped_dict.get(sample_id)
        if ped_sample.father is not None:
            if ped_sample.father not in self.participants:
                self.populate_participants(ped_sample.father)
            if self.ped_dict.get(ped_sample.father).affected:
                affected_parents.add(ped_sample.father)
            else:
                unaffected_parents.add(ped_sample.father)

        if ped_sample.mother is not None:
            if ped_sample.mother not in self.participants:
                self.populate_participants(ped_sample.mother)
            if self.ped_dict.get(ped_sample.mother).affected:
                affected_parents.add(ped_sample.mother)
            else:
                unaffected_parents.add(ped_sample.mother)

        self.participants[sample_id] = Participant(
            details=ped_sample,
            mother=self.participants.get(ped_sample.mother),
            father=self.participants.get(ped_sample.father),
            children=[],
            affected_parents=affected_parents,
            unaffected_parents=unaffected_parents,
        )

    def apply_children(self):
        """
        iterates over the family members and adds children to nodes
        this allows for bi-directional inheritance checks
        :return:
        """

        # flick through all participants
        for participant in self.participants.values():
            # if this person has a father, add child to father
            if participant.father is not None:
                participant.father.children.append(participant)
            # repeat for mother
            if participant.mother is not None:
                participant.mother.children.append(participant)
