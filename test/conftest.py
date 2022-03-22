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


@pytest.fixture(name='cleanup', scope='session', autouse=True)
def fixture_hail_cleanup():
    """
    a fixture to clean up hail log files
    irrelevant in CI, a right pain for local testing
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
