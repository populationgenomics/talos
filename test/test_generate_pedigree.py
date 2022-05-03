"""
methods for testing the sample-metadata API queries
"""

from copy import deepcopy

from unittest.mock import patch

import json
import os

from reanalysis.generate_pedigree import (
    ext_to_int_sample_map,
    get_fat_pedigree,
    get_pedigree_for_project,
)


PWD = os.path.dirname(__file__)
PROJECT = 'fake-project'
INPUT = os.path.join(PWD, 'input')
JSON_PED = os.path.join(INPUT, 'mock_pedigree.json')
LOOKUP_PED = os.path.join(INPUT, 'mock_sm_lookup.json')

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
    },
    {
        'family_id': 'FAMX',
        'individual_id': 'samx',
        'paternal_id': 'samy',
        'maternal_id': 'samz',
        'sex': 1,
        'affected': 1,
    },
]


@patch('reanalysis.generate_pedigree.FamilyApi.get_pedigree')
def test_get_pedigree_for_project(
    get_ped_patch,
):
    """

    :param get_ped_patch:
    :return:
    """

    with open(JSON_PED, 'r', encoding='utf-8') as handle:
        payload = json.load(handle)
        get_ped_patch.return_value = payload
        result = get_pedigree_for_project(project=PROJECT)
        assert result == payload


@patch(
    'reanalysis.generate_pedigree.ParticipantApi.'
    'get_external_participant_id_to_internal_sample_id'
)
def test_ext_to_int_sample_map(
    map_mock,
):
    """
    fetch method using a mocked API endpoint
    :param map_mock:
    :return:
    """

    with open(LOOKUP_PED, 'r', encoding='utf-8') as handle:
        payload = json.load(handle)
        map_mock.return_value = payload
        result = ext_to_int_sample_map(project=PROJECT)
        assert isinstance(result, dict)
        assert result == {
            'FAM1_father': ['CPG11'],
            'FAM1_mother': ['CPG12'],
            'FAM1_proband': ['CPG13'],
            'FAM2_proband': ['CPG41'],
        }


def test_get_clean_pedigree():
    """

    :return:
    """
    cleaned = get_fat_pedigree(
        pedigree_dicts=deepcopy(DIRTY_PED),
        sample_to_cpg_dict=SAMPLE_TO_CPG,
        singles=False,
        plink=False,
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


def test_get_clean_pedigree_singles():
    """

    :return:
    """
    cleaned = get_fat_pedigree(
        pedigree_dicts=deepcopy(DIRTY_PED),
        sample_to_cpg_dict=SAMPLE_TO_CPG,
        singles=True,
        plink=False,
    )
    assert cleaned == [
        {
            'family_id': '1',
            'individual_id': ['cpg1'],
            'paternal_id': [''],
            'maternal_id': [''],
            'sex': 1,
            'affected': 1,
        }
    ]


def test_get_clean_pedigree_singles_plink():
    """

    :return:
    """
    cleaned = get_fat_pedigree(
        pedigree_dicts=deepcopy(DIRTY_PED),
        sample_to_cpg_dict=SAMPLE_TO_CPG,
        singles=True,
        plink=True,
    )
    assert cleaned == [
        {
            'family_id': '1',
            'individual_id': ['cpg1'],
            'paternal_id': ['0'],
            'maternal_id': ['0'],
            'sex': 1,
            'affected': 1,
        }
    ]
