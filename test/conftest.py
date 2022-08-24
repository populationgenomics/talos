"""
A home for common test fixtures
"""

import os
from typing import Any
import pytest

from cyvcf2 import VCFReader
import hail as hl
from hail.utils.java import FatalError
from peddy.peddy import Ped

import cpg_utils
from reanalysis.utils import AbstractVariant, read_json_from_path


PWD = os.path.dirname(__file__)
INPUT = os.path.join(PWD, 'input')

HAIL_VCF = os.path.join(INPUT, 'single_hail.vcf.bgz')
DE_NOVO_TRIO = os.path.join(INPUT, 'de_novo.vcf.bgz')
DE_NOVO_PED = os.path.join(INPUT, 'de_novo_ped.fam')
QUAD_PED = os.path.join(INPUT, 'trio_plus_sibling.fam')
LABELLED = os.path.join(INPUT, '1_labelled_variant.vcf.bgz')
TEST_CONF = os.path.join(INPUT, 'test_conf.json')
PED_FILE = os.path.join(INPUT, 'pedfile.ped')
AIP_OUTPUT = os.path.join(INPUT, 'aip_output_example.json')
SEQR_OUTPUT = os.path.join(INPUT, 'seqr_tags.tsv')
PHASED_TRIO = os.path.join(INPUT, 'phased_trio.vcf.bgz')
PANELAPP_CHANGES = os.path.join(INPUT, 'panel_changes_expected.json')
LOOKUP_PED = os.path.join(INPUT, 'mock_sm_lookup.json')

# can't use a session-scoped monkeypatch, so just set
os.putenv(
    'CPG_CONFIG_PATH',
    os.path.join(os.path.dirname(cpg_utils.__file__), 'config-template.toml'),
)


@pytest.fixture(name='sm_lookup', scope='session')
def fixture_sm_api_lookup() -> dict[str, Any]:
    """
    :return: Ped
    """
    return read_json_from_path(LOOKUP_PED)


@pytest.fixture(name='panel_changes', scope='session')
def fixture_panel_with_new_gene_updates() -> dict[str, Any]:
    """
    :return: Ped
    """
    return read_json_from_path(PANELAPP_CHANGES)


@pytest.fixture(name='peddy_ped', scope='session')
def fixture_peddy_ped() -> Ped:
    """

    :return: Ped
    """
    return Ped(PED_FILE)


@pytest.fixture(name='cleanup', scope='session', autouse=True)
def fixture_hail_cleanup():
    """
    a fixture to clean up hail log files
    irrelevant in CI, a right pain for local testing
    auto-use + session + immediate yield means this is the last method call
    :return:
    """

    # start hail once for the whole runtime
    try:
        hl.init(default_reference='GRCh38')
    except FatalError:
        print('failure - hail already initiated')

    # yield something to suspend
    yield ''

    parent_dir = os.path.join(PWD, os.pardir)

    log_files = [
        filename
        for filename in os.listdir(parent_dir)
        if filename.startswith('hail') and os.path.splitext(filename)[1] == '.log'
    ]

    # remove all hail log files
    for filename in log_files:
        os.remove(os.path.join(parent_dir, filename))


@pytest.fixture(name='hail_matrix', scope='session')
def fixture_hail_matrix():
    """
    loads the single variant as a matrix table
    :return:
    """
    return hl.import_vcf(HAIL_VCF, reference_genome='GRCh38')


@pytest.fixture(name='single_variant_vcf_path')
def fixture_single_variant_vcf_path():
    """
    passes path to the single variant VCF
    :return:
    """

    return HAIL_VCF


@pytest.fixture(name='phased_vcf_path')
def fixture_phased_trio_vcf_path():
    """
    passes path to the phased trio VCF
    :return:
    """

    return PHASED_TRIO


@pytest.fixture(name='phased_variants')
def fixture_phased_trio_variants():
    """
    passes path to the phased trio VCF
    :return:
    """

    conf_json = read_json_from_path(TEST_CONF)
    vcf_reader = VCFReader(PHASED_TRIO)
    two_variants = [
        AbstractVariant(var, vcf_reader.samples, config=conf_json) for var in vcf_reader
    ]
    return two_variants


@pytest.fixture(name='trio_ped')
def fixture_trio_ped():
    """
    sends the location of the Trio Pedigree (PLINK)
    :return:
    """

    return DE_NOVO_PED


@pytest.fixture(name='quad_ped')
def fixture_quad_ped():
    """
    sends the location of the Quad Pedigree (PLINK)
    :return:
    """

    return QUAD_PED


@pytest.fixture(name='conf_json_path')
def fixture_test_conf_path():
    """
    returns the path to the config JSON file
    :return:
    """
    return TEST_CONF


@pytest.fixture(name='trio_abs_variant')
def fixture_trio_abs_variant():
    """
    sends the location of the Trio Pedigree (PLINK)
    Cat. 3, and Cat. 4 for PROBAND only
    :return:
    """
    conf_json = read_json_from_path(TEST_CONF)
    vcf_reader = VCFReader(LABELLED)
    cyvcf_var = next(vcf_reader)

    return AbstractVariant(cyvcf_var, vcf_reader.samples, config=conf_json)


@pytest.fixture(name='two_trio_abs_variants')
def fixture_two_trio_abs_variants():
    """
    sends the location of the Trio Pedigree (PLINK)
    1) Cat. 3, and Cat. 4 for PROBAND only
    2) Cat. 1 + 3, and Cat. 4 for PROBAND only
    :return:
    """
    conf_json = read_json_from_path(TEST_CONF)
    vcf_reader = VCFReader(LABELLED)
    two_variants = [
        AbstractVariant(var, vcf_reader.samples, config=conf_json) for var in vcf_reader
    ]
    return two_variants


@pytest.fixture(name='two_trio_variants_vcf')
def fixture_path_to_two_trio_abs_variants():
    """
    sends the location of the Trio VCF
    :return:
    """
    return LABELLED


@pytest.fixture(name='de_novo_matrix')
def fixture_de_novo_matrix():
    """
    loads the single variant, trio VCF, as a matrix table
    :return:
    """
    return hl.import_vcf(DE_NOVO_TRIO, reference_genome='GRCh38')


@pytest.fixture(name='output_json', scope='session')
def fixture_output_json():
    """
    loads and returns the JSON of the output
    :return:
    """

    return read_json_from_path(AIP_OUTPUT)


@pytest.fixture(name='seqr_csv_output', scope='session')
def fixture_output_seqr_tsv():
    """
    returns the path to the TSV of Seqr variants
    :return:
    """

    return SEQR_OUTPUT


@pytest.fixture(
    params=[
        ('x', 'x', 'frameshift_variant', 'protein_coding', '', 1),
        ('x', 'x', 'frameshift_variant', '', 'NM_relevant', 1),
        ('x', 'o', 'frameshift_variant', 'protein_coding', 'NM_relevant', 0),
        ('x', 'x', 'frameshift_variant', '', '', 0),
    ],
    name='csq_matrix',
    scope='session',
)
def fixture_csq_matrix(request, hail_matrix):
    """

    :param request: the keyword for access to fixture.params
    :param hail_matrix:
    :return:
    """

    gene_ids, gene_id, consequences, biotype, mane_select, row = request.param

    return (
        hail_matrix.annotate_rows(
            geneIds=gene_ids,
            vep=hl.Struct(
                transcript_consequences=hl.array(
                    [
                        hl.Struct(
                            consequence_terms=hl.set([consequences]),
                            biotype=biotype,
                            gene_id=gene_id,
                            mane_select=mane_select,
                        )
                    ]
                ),
            ),
        ),
        row,
    )
