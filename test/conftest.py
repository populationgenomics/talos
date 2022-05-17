"""
A home for common test fixtures
"""


import os
import pytest
import hail as hl
from hail.utils.java import FatalError


PWD = os.path.dirname(__file__)
INPUT = os.path.join(PWD, 'input')

# contains a single variant at chr1:1, with minimal info
HAIL_VCF = os.path.join(INPUT, 'single_hail.vcf.bgz')
HAIL_MULTI_SAM = os.path.join(INPUT, 'multiple_hail.vcf.bgz')
DE_NOVO_TRIO = os.path.join(INPUT, 'de_novo.vcf.bgz')
DE_NOVO_PED = os.path.join(INPUT, 'de_novo_ped.fam')


@pytest.fixture(name='hail_matrix')
def fixture_hail_matrix():
    """
    loads the single variant as a matrix table
    :return:
    """
    try:
        hl.init(default_reference='GRCh38')
    except FatalError:
        print('failure - hail already initiated')
    return hl.import_vcf(HAIL_VCF, reference_genome='GRCh38')


@pytest.fixture(name='trio_ped')
def fixture_trio_ped():
    """
    sends the location of the Trio Pedigree (PLINK)
    :return:
    """

    return DE_NOVO_PED


@pytest.fixture(name='de_novo_matrix')
def fixture_de_novo_matrix():
    """
    loads the single variant, trio VCF, as a matrix table
    :return:
    """
    try:
        hl.init(default_reference='GRCh38')
    except FatalError:
        print('failure - hail already initiated')
    return hl.import_vcf(DE_NOVO_TRIO, reference_genome='GRCh38')


@pytest.fixture(name='hail_comp_het')
def fixture_hail_matrix_comp_het():
    """
    loads the single variant as a matrix table
    :return:
    """
    try:
        hl.init(default_reference='GRCh38')
    except FatalError:
        print('failure - hail already initiated')
    return hl.import_vcf(HAIL_MULTI_SAM, reference_genome='GRCh38')


@pytest.fixture(name='cleanup', scope='session', autouse=True)
def fixture_hail_cleanup():
    """
    a fixture to clean up hail log files
    irrelevant in CI, a right pain for local testing
    auto-use + session + immediate yield means this is the last method call
    :return:
    """
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
