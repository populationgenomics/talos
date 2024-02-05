"""
methods for testing the sample-metadata API queries
"""

from copy import deepcopy
from unittest.mock import patch

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

    Args:
        map_mock ():
        gql_mock ():
        sm_lookup ():
    """
    gql_mock.return_value = None
    map_mock.return_value = sm_lookup
    result = ext_to_int_sample_map(project=PROJECT)
    assert isinstance(result, dict)
    assert result == {
        'FAM1_father': ['CPT11'],
        'FAM1_mother': ['CPT12'],
        'FAM1_proband': ['CPT13'],
        'FAM2_proband': ['CPT41'],
    }


def test_get_clean_pedigree_fails(caplog):
    """
    this no longer fails, but does report an error
    """
    ped = [{'individual_id': 'plink'}]

    dict_list = get_ped_with_permutations(pedigree_dicts=ped, ext_lookup={})
    assert not dict_list
    assert 'plink' in caplog.text


def test_get_clean_pedigree():
    """

    :return:
    """
    cleaned = get_ped_with_permutations(
        pedigree_dicts=deepcopy(DIRTY_PED), ext_lookup=SAMPLE_TO_CPG
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
