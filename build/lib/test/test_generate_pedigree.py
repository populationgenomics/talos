"""
methods for testing the sample-metadata API queries
"""

from copy import deepcopy
from unittest.mock import patch
import pytest

from helpers.prepare_aip_cohort import get_ped_with_permutations
from helpers.utils import ext_to_int_sample_map

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


@patch('helpers.utils.gql')
@patch('helpers.utils.query')
def test_ext_to_int_sample_map(map_mock, gql_mock, sm_lookup):
    """
    fetch method using a mocked API endpoint
    :param map_mock:
    :return:
    """
    gql_mock.return_value = None
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

    with pytest.raises(KeyError):
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
