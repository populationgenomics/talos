"""
A home for common test fixtures
"""

from os import environ
from os.path import join
from pathlib import Path
from typing import Any

import pytest
from _pytest.logging import LogCaptureFixture
from cyvcf2 import VCFReader
from loguru import logger

import hail as hl

from talos.data_model import BaseFields, Entry, SneakyTable, TXFields, VepVariant
from talos.utils import create_small_variant, read_json_from_path

# force this to come first
PWD = Path(__file__).parent
INPUT: str = str(PWD / 'input')
hl.init()
hl.default_reference('GRCh38')
environ['TALOS_CONFIG'] = join(INPUT, 'config.toml')

LABELLED = join(INPUT, '1_labelled_variant.vcf.bgz')
Talos_OUTPUT = join(INPUT, 'aip_output_example.json')
DE_NOVO_PED = join(INPUT, 'de_novo_ped.fam')
FAKE_OBO = join(INPUT, 'hpo_test.obo')
LOOKUP_PED = join(INPUT, 'mock_sm_lookup.json')
PHASED_TRIO = join(INPUT, 'newphase.vcf.bgz')
PED_FILE = join(INPUT, 'pedfile.ped')
SEQR_OUTPUT = join(INPUT, 'seqr_tags.tsv')
QUAD_PED = join(INPUT, 'trio_plus_sibling.fam')
SUB_STUB = join(INPUT, 'tiny_summary.txt.gz')

# panelapp testing paths
PANEL_ACTIVITIES = join(INPUT, 'panelapp_activities.json')
PANELAPP_LATEST = join(INPUT, 'panelapp_current_137.json')
PANELAPP_ALL_PANELS = join(INPUT, 'panelapp_all_panels.json')
PANELAPP_INCIDENTALOME = join(INPUT, 'incidentalome.json')
FAKE_PANELAPP_OVERVIEW = join(INPUT, 'panel_overview.json')


@pytest.fixture
def caplog(caplog: LogCaptureFixture):
    handler_id = logger.add(
        caplog.handler,
        format='{message}',
        level=0,
        filter=lambda record: record['level'].no >= caplog.handler.level,
        enqueue=False,  # Set to 'True' if your test is spawning child processes.
    )
    yield caplog
    logger.remove(handler_id)


@pytest.fixture(name='cleanup', scope='session', autouse=True)
def fixture_hail_cleanup():
    """
    a fixture to clean up hail log files
    irrelevant in CI, a right pain for local testing
    auto-use + session + immediate yield means this is the last method call
    """

    # yield something to suspend
    yield ''

    log_files = [
        filename for filename in PWD.parent.iterdir() if filename.name.startswith('hail') and filename.suffix == '.log'
    ]

    # remove all hail log files
    for filename in log_files:
        filename.unlink()


@pytest.fixture(name='make_a_mt', scope='session')
def fixture_make_a_mt(tmp_path_factory) -> hl.MatrixTable:
    """
    a fixture to make a matrix table
    """
    tmp_path = tmp_path_factory.mktemp('mt_goes_here')
    sample_gt = Entry('0/1')
    sample_data = {'SAMPLE': sample_gt}
    sample_schema = {'SAMPLE': sample_gt.get_schema_entry()}
    v = VepVariant(
        BaseFields('chr1:12345', alleles=['A', 'G']),
        [TXFields('a', 'ensga')],
        sample_data=sample_data,
    )
    return SneakyTable([v], sample_details=sample_schema, tmp_path=str(tmp_path)).to_hail()


@pytest.fixture(name='make_a_vcf', scope='session')
def fixture_make_a_vcf(make_a_mt, tmp_path_factory) -> str:
    """
    a fixture to make a matrix table
    """
    tmp_path = tmp_path_factory.mktemp('vcf_goes_here')
    vcf_path = str(tmp_path / 'test.vcf.bgz')
    hl.export_vcf(make_a_mt, vcf_path, tabix=True)
    return vcf_path


@pytest.fixture(name='test_input_path', scope='session')
def fixture_test_input_path() -> str:
    """path to the test input directory"""
    return INPUT


@pytest.fixture(name='test_input_models_path', scope='session')
def fixture_test_input_models_path() -> str:
    """path to the test input directory"""
    return join(INPUT, 'models')


@pytest.fixture(name='fake_obo_path', scope='session')
def fixture_fake_obo() -> str:
    """path to fake obo"""
    return FAKE_OBO


@pytest.fixture(name='fake_panelapp_overview', scope='session')
def fixture_panelapp_overview() -> Any:
    """path to panelapp_overview json"""
    return read_json_from_path(FAKE_PANELAPP_OVERVIEW)


@pytest.fixture(name='panel_activities', scope='session')
def fixture_panel_activities() -> Any:
    """path to activities json"""
    return read_json_from_path(PANEL_ACTIVITIES)


@pytest.fixture(name='latest_mendeliome', scope='session')
def fixture_latest_mendeliome() -> Any:
    """path to incidentalome json"""
    return read_json_from_path(PANELAPP_LATEST)


@pytest.fixture(name='panels_and_hpos', scope='session')
def fixture_panels_and_hpo() -> Any:
    return read_json_from_path(PANELAPP_ALL_PANELS)


@pytest.fixture(name='latest_incidentalome', scope='session')
def fixture_latest_incidentalome() -> Any:
    """path to mendeliome json"""
    return read_json_from_path(PANELAPP_INCIDENTALOME)


@pytest.fixture(name='sm_lookup', scope='session')
def fixture_sm_api_lookup() -> Any:
    """
    :return: Ped
    """
    return read_json_from_path(LOOKUP_PED)


@pytest.fixture(name='pedigree_path', scope='session')
def pedigree_path() -> str:
    """
    :return: Ped
    """
    return PED_FILE


@pytest.fixture(name='phased_vcf_path')
def fixture_phased_trio_vcf_path():
    """path to the phased trio VCF"""

    return PHASED_TRIO


@pytest.fixture(name='phased_variants')
def fixture_phased_trio_variants():
    """path to the phased trio VCF"""

    vcf_reader = VCFReader(PHASED_TRIO)
    return [create_small_variant(var, vcf_reader.samples) for var in vcf_reader]


@pytest.fixture(name='trio_ped')
def fixture_trio_ped():
    """location of the Trio Pedigree (PLINK)"""

    return DE_NOVO_PED


@pytest.fixture(name='quad_ped')
def fixture_quad_ped():
    """location of the Quad Pedigree (PLINK)"""

    return QUAD_PED


@pytest.fixture(name='trio_abs_variant')
def fixture_trio_abs_variant():
    """
    sends the location of the Trio Pedigree (PLINK)
    Cat. 3, and Cat. 4 for PROBAND only
    """
    vcf_reader = VCFReader(LABELLED)
    cyvcf_var = next(vcf_reader)

    return create_small_variant(cyvcf_var, vcf_reader.samples)


@pytest.fixture(name='cyvcf_example_variant')
def fixture_cyvcf_variant():
    """
    sends the location of the Trio Pedigree (PLINK)
    Cat. 3, and Cat. 4 for PROBAND only
    """
    vcf_reader = VCFReader(LABELLED)
    return next(vcf_reader)


@pytest.fixture(name='two_trio_abs_variants')
def fixture_two_trio_abs_variants():
    """
    sends the location of the Trio Pedigree (PLINK)
    1) Cat. 3, and Cat. 4 for PROBAND only
    2) Cat. 1 + 3, and Cat. 4 for PROBAND only
    """
    vcf_reader = VCFReader(LABELLED)
    return [create_small_variant(var, vcf_reader.samples) for var in vcf_reader]


@pytest.fixture(name='two_trio_variants_vcf')
def fixture_path_to_two_trio_abs_variants():
    """sends the location of the Trio VCF"""
    return LABELLED


@pytest.fixture(name='output_json', scope='session')
def fixture_output_json():
    """returns dict of the JSON output"""

    return read_json_from_path(Talos_OUTPUT)


@pytest.fixture(name='seqr_csv_output', scope='session')
def fixture_output_seqr_tsv():
    """path to the TSV of Seqr variants"""

    return SEQR_OUTPUT


@pytest.fixture(name='sub_stub', scope='session')
def fixture_sub_stub():
    """path to the TXT file of submissions"""

    return SUB_STUB
