"""
A home for common test fixtures
"""

from typing import Any
import pytest

from cyvcf2 import VCFReader
import hail as hl
from peddy.peddy import Ped

from cpg_utils import to_path
from cpg_utils.config import set_config_paths

from reanalysis.data_model import BaseFields, Entry, TXFields, VepVariant, SneakyTable

PWD = to_path(__file__).parent
INPUT = PWD / 'input'

# force this to come first
CONF_BASE = INPUT / 'reanalysis_global.toml'
CONF_COHORT = INPUT / 'reanalysis_cohort.toml'

hl.init(default_reference='GRCh38')
set_config_paths([str(CONF_BASE), str(CONF_COHORT)])

# pylint: disable=wrong-import-position
from reanalysis.utils import AbstractVariant, read_json_from_path

LABELLED = INPUT / '1_labelled_variant.vcf.bgz'
AIP_OUTPUT = INPUT / 'aip_output_example.json'
DE_NOVO_PED = INPUT / 'de_novo_ped.fam'
FAKE_OBO = INPUT / 'hpo_test.obo'
LOOKUP_PED = INPUT / 'mock_sm_lookup.json'
PHASED_TRIO = INPUT / 'newphase.vcf.bgz'
PED_FILE = INPUT / 'pedfile.ped'
SEQR_OUTPUT = INPUT / 'seqr_tags.tsv'
QUAD_PED = INPUT / 'trio_plus_sibling.fam'
SUB_STUB = INPUT / 'tiny_summary.txt.gz'

# panelapp testing paths
PANELAPP_LATEST = INPUT / 'panelapp_current_137.json'
PANELAPP_INCIDENTALOME = INPUT / 'incidentalome.json'
FAKE_PANELAPP_OVERVIEW = INPUT / 'panel_overview.json'


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
        filename
        for filename in PWD.parent.iterdir()
        if filename.name.startswith('hail') and filename.suffix == '.log'
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
    return SneakyTable(
        [v], sample_details=sample_schema, tmp_path=str(tmp_path)
    ).to_hail()


@pytest.fixture(name='make_a_vcf', scope='session')
def fixture_make_a_vcf(make_a_mt, tmp_path_factory) -> str:
    """
    a fixture to make a matrix table
    """
    tmp_path = tmp_path_factory.mktemp('vcf_goes_here')
    vcf_path = str(tmp_path / 'test.vcf.bgz')
    hl.export_vcf(make_a_mt, vcf_path, tabix=True)
    return vcf_path


@pytest.fixture(name='fake_obo_path', scope='session')
def fixture_fake_obo() -> str:
    """path to fake obo"""
    return FAKE_OBO


@pytest.fixture(name='fake_panelapp_overview', scope='session')
def fixture_panelapp_overview() -> Any:
    """path to panelapp_overview json"""
    return read_json_from_path(FAKE_PANELAPP_OVERVIEW)


@pytest.fixture(name='latest_mendeliome', scope='session')
def fixture_latest_mendeliome() -> Any:
    """path to incidentalome json"""
    return read_json_from_path(PANELAPP_LATEST)


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


@pytest.fixture(name='peddy_ped', scope='session')
def fixture_peddy_ped() -> Ped:
    """
    :return: Ped
    """
    return Ped(str(PED_FILE))


@pytest.fixture(name='phased_vcf_path')
def fixture_phased_trio_vcf_path():
    """path to the phased trio VCF"""

    return PHASED_TRIO


@pytest.fixture(name='phased_variants')
def fixture_phased_trio_variants():
    """path to the phased trio VCF"""

    vcf_reader = VCFReader(PHASED_TRIO)
    two_variants = [AbstractVariant(var, vcf_reader.samples) for var in vcf_reader]
    return two_variants


@pytest.fixture(name='trio_ped')
def fixture_trio_ped():
    """location of the Trio Pedigree (PLINK)"""

    return DE_NOVO_PED


@pytest.fixture(name='quad_ped')
def fixture_quad_ped():
    """location of the Quad Pedigree (PLINK)"""

    return Ped(str(QUAD_PED))


@pytest.fixture(name='trio_abs_variant')
def fixture_trio_abs_variant():
    """
    sends the location of the Trio Pedigree (PLINK)
    Cat. 3, and Cat. 4 for PROBAND only
    """
    vcf_reader = VCFReader(LABELLED)
    cyvcf_var = next(vcf_reader)

    return AbstractVariant(cyvcf_var, vcf_reader.samples)


@pytest.fixture(name='two_trio_abs_variants')
def fixture_two_trio_abs_variants():
    """
    sends the location of the Trio Pedigree (PLINK)
    1) Cat. 3, and Cat. 4 for PROBAND only
    2) Cat. 1 + 3, and Cat. 4 for PROBAND only
    """
    vcf_reader = VCFReader(LABELLED)
    two_variants = [AbstractVariant(var, vcf_reader.samples) for var in vcf_reader]
    return two_variants


@pytest.fixture(name='two_trio_variants_vcf')
def fixture_path_to_two_trio_abs_variants():
    """sends the location of the Trio VCF"""
    return LABELLED


@pytest.fixture(name='output_json', scope='session')
def fixture_output_json():
    """returns dict of the JSON output"""

    return read_json_from_path(AIP_OUTPUT)


@pytest.fixture(name='seqr_csv_output', scope='session')
def fixture_output_seqr_tsv():
    """path to the TSV of Seqr variants"""

    return SEQR_OUTPUT


@pytest.fixture(name='sub_stub', scope='session')
def fixture_sub_stub():
    """path to the TXT file of submissions"""

    return SUB_STUB


# @pytest.fixture(scope='session')
# def session_temp_dir(tmp_path_factory):
#     """
#     session-scoped template path
#     return a tempdir base
#     disparate tests can share a tempdir
#     """
#
#     return tmp_path_factory.mktemp('TEMP')
