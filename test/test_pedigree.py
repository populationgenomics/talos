"""
test methods for linked-list pedigree parser
"""

import os

from reanalysis.pedigree import PedigreeParser


PWD = os.path.dirname(__file__)
INPUT = os.path.join(PWD, 'input')
PED_FILE = os.path.join(INPUT, 'pedfile.ped')

# test pedigree consists of a male family (1) and female family (2)


def test_ped_parsing():
    """
    read pedigree from a test file

    :return:
    """

    expected_samples = [
        'father_1',
        'father_2',
        'female',
        'male',
        'mother_1',
        'mother_2',
    ]

    parsed = PedigreeParser(PED_FILE)

    assert sorted(list(parsed.participants.keys())) == expected_samples

    assert not parsed.participants['male'].details.is_female
    assert parsed.participants['male'].details.affected
    assert parsed.participants['male'].father.details.sample_id == 'father_1'
    assert parsed.participants['male'].mother.details.sample_id == 'mother_1'
    assert parsed.participants['mother_1'].children[0].details.sample_id == 'male'

    assert parsed.participants['female'].details.is_female
    assert parsed.participants['female'].details.affected
    assert parsed.participants['female'].children == []
    assert parsed.participants['female'].father.details.sample_id == 'father_2'
    assert parsed.participants['female'].mother.details.sample_id == 'mother_2'
    assert parsed.participants['father_2'].children[0].details.sample_id == 'female'
