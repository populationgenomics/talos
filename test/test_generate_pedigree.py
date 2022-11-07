"""
methods for testing the sample-metadata API queries
"""

from copy import deepcopy
from unittest.mock import patch
import pytest

from helpers.prepare_aip_cohort import (
    ext_to_int_sample_map,
    get_ped_with_permutations,
    hash_reduce_dicts,
)

PROJECT = 'fake-project'

SAMPLE_TO_CPG = {
    'sam1': ['cpg1'],
    'sam2': ['cpg2'],
    'sam3': ['cpg3'],
}
DIRTY_PED = [
    {
        'family_id': 'FAM1',
        'individual_id': 'sam1',
        'paternal_id': 'sam2',
        'maternal_id': 'sam3',
        'sex': 1,
        'affected': 1,
    }
]


@patch(
    'helpers.prepare_aip_cohort.ParticipantApi.'
    'get_external_participant_id_to_internal_sample_id'
)
def test_ext_to_int_sample_map(map_mock, sm_lookup):
    """
    fetch method using a mocked API endpoint
    :param map_mock:
    :return:
    """

    map_mock.return_value = sm_lookup
    result = ext_to_int_sample_map(project=PROJECT)
    assert isinstance(result, dict)
    assert result == {
        'FAM1_father': ['CPG11'],
        'FAM1_mother': ['CPG12'],
        'FAM1_proband': ['CPG13'],
        'FAM2_proband': ['CPG41'],
    }


def test_get_clean_pedigree_fails():
    """

    :return:
    """
    ped = [{'individual_id': 'plink'}]

    with pytest.raises(Exception):
        get_ped_with_permutations(pedigree_dicts=ped, sample_to_cpg_dict={})


def test_get_clean_pedigree():
    """

    :return:
    """
    cleaned = get_ped_with_permutations(
        pedigree_dicts=deepcopy(DIRTY_PED), sample_to_cpg_dict=SAMPLE_TO_CPG
    )
    assert cleaned == [
        {
            'family_id': 'FAM1',
            'individual_id': ['cpg1'],
            'paternal_id': ['cpg2'],
            'maternal_id': ['cpg3'],
            'sex': 1,
            'affected': 1,
        }
    ]


def test_hash_reduce():
    """
    examples of the int hash % 100
    family0 35
    family1 76
    family2 97
    family3 76
    family4 55
    family5 43
    family6 35

    note collision of family1 & 3, 0 & 6; not a problem
    :return:
    """
    # test input includes duplication
    pedigree_dicts = [
        {'family_id': 'family0'},
        {'family_id': 'family1'},
        {'family_id': 'family1'},
        {'family_id': 'family2'},
        {'family_id': 'family3'},
        {'family_id': 'family4'},
        {'family_id': 'family4'},
    ]
    result_1 = hash_reduce_dicts(pedigree_dicts, 76)
    assert len(result_1) == 3
    fam1 = [x for x in result_1 if x['family_id'] == 'family1']
    assert len(fam1) == 0

    result_2 = hash_reduce_dicts(pedigree_dicts, 77)
    assert len(result_2) == 6

    # check all of family1 were retained
    fam1 = [x for x in result_2 if x['family_id'] == 'family1']
    assert len(fam1) == 2
