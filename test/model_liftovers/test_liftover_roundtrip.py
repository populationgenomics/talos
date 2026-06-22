"""
Round-trip regression tests for the model liftover engine.

These lock in the *current* behaviour of ``lift_up_model_version`` across every
historical model version, so that the planned migration-engine refactor (see
MODEL_MIGRATION_PLAN.md, Track A) can be shown to be behaviour-preserving.

For ``ResultData`` (the longest/most complex chain) we assert that:
  1. the chain walks up to ``CURRENT_VERSION``
  2. the lifted dict validates against the model
  3. the lifted dict matches a committed golden snapshot

The lifted ``ResultData`` dict is deterministic, so a golden snapshot is the
strongest regression net. ``DownloadedPanelApp`` / ``PanelApp`` chains inject the
current date during liftover, so those are covered with explicit value
assertions instead of frozen snapshots.

Golden snapshots are generated on first run (when absent) and then committed;
thereafter the test compares against them and fails on any drift.
"""

import json
from glob import glob
from os import makedirs
from os.path import basename, dirname, exists, join

import pytest

from talos.models import (
    CURRENT_VERSION,
    Coordinates,
    DownloadedPanelApp,
    PanelApp,
    ReportVariant,
    ResultData,
    SmallVariant,
    StructuralVariant,
    lift_up_model_version,
)

GOLDEN_SUBDIR = 'golden'


def _result_data_fixtures(models_path: str) -> list[str]:
    """every committed historical ResultData fixture"""
    return sorted(glob(join(models_path, 'result_data_v*.json')))


def _assert_matches_golden(models_path: str, name: str, lifted: dict) -> None:
    """
    Compare lifted output to a committed golden snapshot.

    When the snapshot is absent it is generated and the run passes, so the first
    invocation produces a file to inspect and commit. Subsequent runs compare
    against the committed snapshot and fail on any drift.
    """
    path = join(models_path, GOLDEN_SUBDIR, name)
    serialised = json.dumps(lifted, indent=2, sort_keys=True, default=str) + '\n'

    if not exists(path):
        makedirs(dirname(path), exist_ok=True)
        with open(path, 'w', encoding='utf-8') as handle:
            handle.write(serialised)
        return

    with open(path, encoding='utf-8') as handle:
        expected = handle.read()
    assert serialised == expected, f'lifted output drifted from golden snapshot {name}'


@pytest.fixture(name='result_data_fixture_ids')
def fixture_result_data_fixture_ids(test_input_models_path) -> list[str]:
    """basenames of all ResultData fixtures - used to label parametrised cases"""
    return [basename(p) for p in _result_data_fixtures(test_input_models_path)]


def test_result_data_fixtures_exist(test_input_models_path):
    """guard against the fixture glob silently matching nothing"""
    assert _result_data_fixtures(test_input_models_path), 'no ResultData fixtures discovered'


def _load_result_data_cases(models_path: str) -> list[tuple[str, dict]]:
    cases = []
    for path in _result_data_fixtures(models_path):
        with open(path, encoding='utf-8') as handle:
            cases.append((basename(path), json.load(handle)))
    return cases


def test_result_data_roundtrip(test_input_models_path):
    """
    every historical ResultData fixture must lift to CURRENT_VERSION, validate,
    and match its golden snapshot
    """
    for name, data in _load_result_data_cases(test_input_models_path):
        lifted = lift_up_model_version(data, model=ResultData)

        assert lifted['version'] == CURRENT_VERSION, f'{name} did not lift to {CURRENT_VERSION}'

        parsed = ResultData.model_validate(lifted)
        assert parsed.version == CURRENT_VERSION

        golden_name = basename(name).replace('.json', '.lifted.json')
        _assert_matches_golden(test_input_models_path, golden_name, lifted)


def test_downloaded_panelapp_roundtrip_from_2_1_0():
    """
    a 2.1.0 DownloadedPanelApp must lift through 2.2.0 (date injection) and 2.3.0
    (confidence backfill) and validate
    """
    data = {
        'version': '2.1.0',
        'genes': {
            'ENSG001': {
                'symbol': 'GENE1',
                'chrom': '1',
                'mane_symbol': '',
                'ensg': 'ENSG001',
                'panels': {
                    137: {'moi': 'biallelic', 'date': '2023-01-01'},
                },
            },
        },
        'versions': [],
        'hpos': {},
    }
    lifted = lift_up_model_version(data, model=DownloadedPanelApp)

    assert lifted['version'] == CURRENT_VERSION
    # 2.1.0 -> 2.2.0 backfills a download date
    assert lifted['date']
    parsed = DownloadedPanelApp.model_validate(lifted)
    # 2.2.0 -> 2.3.0 backfills per-panel confidence of 3
    assert parsed.genes['ENSG001'].panels[137].confidence == 3


def test_panelapp_1_2_0_to_2_0_0_is_not_liftable():
    """
    PanelApp is fundamentally different across the 1.2.0 -> 2.0.0 boundary; the
    engine must refuse rather than silently produce a broken object
    """
    data = {'version': '1.2.0', 'metadata': {}, 'genes': {}, 'participants': {}}
    with pytest.raises(ValueError, match='regenerate'):
        lift_up_model_version(data, model=PanelApp)


def test_variant_union_discriminates_on_kind():
    """
    a ReportVariant must deserialise its var_data back into the correct subtype via the
    `kind` discriminator - the StructuralVariant case is not covered by any liftover fixture
    """
    coords = Coordinates(chrom='1', pos=1, ref='A', alt='C')

    for builder, expected in (
        (SmallVariant(coordinates=coords, transcript_consequences=[]), SmallVariant),
        (StructuralVariant(coordinates=coords), StructuralVariant),
    ):
        report = ReportVariant(sample='sam1', var_data=builder)
        reloaded = ReportVariant.model_validate(report.model_dump())
        assert isinstance(reloaded.var_data, expected)


def test_panelapp_roundtrip_from_2_0_0():
    """a 2.0.0 PanelApp must drop the per-participant ext_id, lift, and validate"""
    data = {
        'version': '2.0.0',
        'metadata': {137: {'id': 137, 'name': 'panel', 'version': '1'}},
        'genes': {
            'GENE1': {'symbol': 'GENE1', 'chrom': '1', 'moi': 'Biallelic', 'new': [], 'panels': [137]},
        },
        'participants': {
            'sam1': {
                'ext_id': 'EXT1',
                'family_id': '1',
                'hpo_terms': [],
                'panels': [137],
                'matched_genes': [],
                'matched_phenotypes': [],
            },
        },
    }
    lifted = lift_up_model_version(data, model=PanelApp)

    assert lifted['version'] == CURRENT_VERSION
    assert 'ext_id' not in lifted['participants']['sam1']
    parsed = PanelApp.model_validate(lifted)
    assert 'sam1' in parsed.participants
