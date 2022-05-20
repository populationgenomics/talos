"""
A home for common test fixtures
"""


import json
import os
import pytest
import hail as hl
from hail.utils.java import FatalError

from cyvcf2 import VCFReader
from peddy.peddy import Ped

from reanalysis.utils import AbstractVariant


PWD = os.path.dirname(__file__)
INPUT = os.path.join(PWD, 'input')

# contains a single variant at chr1:1, with minimal info
HAIL_VCF = os.path.join(INPUT, 'single_hail.vcf.bgz')
HAIL_MULTI_SAM = os.path.join(INPUT, 'multiple_hail.vcf.bgz')
DE_NOVO_TRIO = os.path.join(INPUT, 'de_novo.vcf.bgz')
DE_NOVO_PED = os.path.join(INPUT, 'de_novo_ped.fam')
LABELLED = os.path.join(INPUT, '1_labelled_variant.vcf.bgz')
TEST_CONF = os.path.join(INPUT, 'test_conf.json')
PED_FILE = os.path.join(INPUT, 'pedfile.ped')


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


@pytest.fixture(name='hail_matrix')
def fixture_hail_matrix():
    """
    loads the single variant as a matrix table
    :return:
    """
    return hl.import_vcf(HAIL_VCF, reference_genome='GRCh38')


@pytest.fixture(name='trio_ped')
def fixture_trio_ped():
    """
    sends the location of the Trio Pedigree (PLINK)
    :return:
    """

    return DE_NOVO_PED


@pytest.fixture(name='trio_abs_variant')
def fixture_trio_abs_variant():
    """
    sends the location of the Trio Pedigree (PLINK)
    Cat. 3, and Cat. 4 for PROBAND only
    :return:
    """

    with open(TEST_CONF, 'r', encoding='utf-8') as handle:
        conf_json = json.load(handle)

    vcf_reader = VCFReader(LABELLED)
    cyvcf_var = next(vcf_reader)

    return AbstractVariant(cyvcf_var, vcf_reader.samples, config=conf_json)


@pytest.fixture(name='de_novo_matrix')
def fixture_de_novo_matrix():
    """
    loads the single variant, trio VCF, as a matrix table
    :return:
    """
    return hl.import_vcf(DE_NOVO_TRIO, reference_genome='GRCh38')


@pytest.fixture(name='hail_comp_het')
def fixture_hail_matrix_comp_het():
    """
    loads the single variant as a matrix table
    :return:
    """
    return hl.import_vcf(HAIL_MULTI_SAM, reference_genome='GRCh38')
