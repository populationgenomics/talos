"""
Unit tests for generate_dataset_prefix in cpg_flow_utils.

The function depends on the CPG config singleton (config_retrieve / dataset_path)
and on cpg_flow Cohort objects, so all external touch-points are mocked. The goal
is full line coverage of generate_dataset_prefix.
"""

from unittest.mock import patch

import pytest

from talos.cpg_internal_scripts import cpg_flow_utils
from talos.cpg_internal_scripts.cpg_flow_utils import generate_dataset_prefix


class FakeDataset:
    """minimal stand-in for a cpg_flow Dataset (only .name is read)"""

    def __init__(self, name: str):
        self.name = name


class FakeCohort:
    """
    minimal stand-in for a cpg_flow Cohort exposing .id and .dataset.name.
    Must be hashable (default object identity hash) because generate_dataset_prefix
    is @cache decorated and receives the cohort as a key.
    """

    def __init__(self, cohort_id: str, dataset_name: str):
        self.id = cohort_id
        self.dataset = FakeDataset(dataset_name)


def make_cohort(cohort_id: str = 'COH123', dataset_name: str = 'cohort-dataset'):
    """build a stand-in for a cpg_flow Cohort with the attributes the function reads"""
    return FakeCohort(cohort_id, dataset_name)


def config_side_effect(seq_type: str = 'genome', long_read: bool = False):
    """
    return a side_effect callable mimicking config.config_retrieve for the two
    keys the function queries
    """

    def _retrieve(keys, default=None):
        if keys == ['workflow', 'sequencing_type']:
            return seq_type
        if keys == ['workflow', 'long_read']:
            return long_read
        raise AssertionError(f'unexpected config key requested: {keys}')

    return _retrieve


@pytest.fixture(autouse=True)
def _clear_cache():
    """generate_dataset_prefix is @cache decorated - reset between tests"""
    generate_dataset_prefix.cache_clear()
    yield
    generate_dataset_prefix.cache_clear()


def test_raises_when_no_cohort_or_dataset():
    """both cohort and dataset None -> RuntimeError"""
    with pytest.raises(RuntimeError, match='Must populate either cohort or dataset'):
        generate_dataset_prefix()


def test_cohort_path_standard():
    """
    cohort supplied (short-read, genome): cohort_id and dataset are pulled from the
    Cohort, no exome/long_read elements, just 'talos/<cohort_id>'
    """
    cohort = make_cohort(cohort_id='COHORTX', dataset_name='ds-from-cohort')

    with (
        patch.object(cpg_flow_utils.config, 'config_retrieve', side_effect=config_side_effect()),
        patch.object(cpg_flow_utils.config, 'dataset_path', return_value='gs://bucket/talos/COHORTX') as ds_path,
        patch.object(cpg_flow_utils, 'to_path', side_effect=lambda x: f'PATH::{x}') as to_path,
    ):
        result = generate_dataset_prefix(cohort=cohort)

    ds_path.assert_called_once_with(suffix='talos/COHORTX', dataset='ds-from-cohort', category=None)
    to_path.assert_called_once_with('gs://bucket/talos/COHORTX')
    assert result == 'PATH::gs://bucket/talos/COHORTX'


def test_dataset_only_no_cohort_id():
    """
    dataset supplied without a cohort (index-page case): cohort_id is None and is
    therefore dropped from the suffix
    """
    with (
        patch.object(cpg_flow_utils.config, 'config_retrieve', side_effect=config_side_effect()),
        patch.object(cpg_flow_utils.config, 'dataset_path', return_value='gs://b/talos') as ds_path,
        patch.object(cpg_flow_utils, 'to_path', side_effect=lambda x: x),
    ):
        result = generate_dataset_prefix(dataset='plain-dataset', category='analysis')

    ds_path.assert_called_once_with(suffix='talos', dataset='plain-dataset', category='analysis')
    assert result == 'gs://b/talos'


def test_exome_and_long_read_and_hash_and_stage():
    """
    exercise the exome and long_read branches plus hash_value and stage_name so every
    element of the ordered suffix list is a non-None string
    """
    cohort = make_cohort(cohort_id='CID', dataset_name='ds')

    with (
        patch.object(cpg_flow_utils.config, 'config_retrieve', side_effect=config_side_effect('exome', True)),
        patch.object(cpg_flow_utils.config, 'dataset_path', return_value='out') as ds_path,
        patch.object(cpg_flow_utils, 'to_path', side_effect=lambda x: x),
    ):
        generate_dataset_prefix(cohort=cohort, hash_value='abc123', stage_name='MyStage')

    ds_path.assert_called_once_with(
        suffix='long_read/exome/talos/CID/abc123/MyStage',
        dataset='ds',
        category=None,
    )


def test_long_read_only():
    """long_read set but genome sequencing -> long_read element present, no exome"""
    with (
        patch.object(cpg_flow_utils.config, 'config_retrieve', side_effect=config_side_effect('genome', True)),
        patch.object(cpg_flow_utils.config, 'dataset_path', return_value='out') as ds_path,
        patch.object(cpg_flow_utils, 'to_path', side_effect=lambda x: x),
    ):
        generate_dataset_prefix(dataset='ds', stage_name='S')

    ds_path.assert_called_once_with(suffix='long_read/talos/S', dataset='ds', category=None)


def test_exome_only():
    """exome sequencing but short read -> exome element present, no long_read"""
    with (
        patch.object(cpg_flow_utils.config, 'config_retrieve', side_effect=config_side_effect('exome', False)),
        patch.object(cpg_flow_utils.config, 'dataset_path', return_value='out') as ds_path,
        patch.object(cpg_flow_utils, 'to_path', side_effect=lambda x: x),
    ):
        generate_dataset_prefix(dataset='ds', hash_value='HSH')

    ds_path.assert_called_once_with(suffix='exome/talos/HSH', dataset='ds', category=None)
