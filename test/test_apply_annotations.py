"""
tests for the annotation module

pytest in GitHub can't read the GS paths - needs local data
or mocking?
"""


import string
import time
import os
from random import choices
import hail as hl
import pytest
from cpg_utils.hail_batch import dataset_path, reference_path

from reanalysis.annotation import apply_annotations


PWD = os.path.dirname(__file__)
INPUT = os.path.join(PWD, 'input')


def timestamp():
    """
    Generate a timestamp plus a short random string for guaranteed uniqueness.
    """
    rand_bit = ''.join(choices(string.ascii_uppercase + string.digits, k=3))
    return time.strftime('%Y-%m%d-%H%M') + rand_bit


INTERVAL = 'chr20-5111495-5111607'


@pytest.mark.gcp
def test_apply_annotations():
    """
    Test annotation.apply_annotations: convert VCF to a MatrixTable,
    and add VEP and other annotations.
    """
    vcf_path = dataset_path(f'unittest/inputs/chr20/joint-called-{INTERVAL}.vcf.gz')
    vep_ht_path = dataset_path(f'unittest/inputs/chr20/vep/chr20-5111495-5111607.ht')
    tmp_bucket = dataset_path(
        f'unittest/test_apply_annotations/{timestamp()}', category='tmp'
    )
    out_mt_path = f'{tmp_bucket}/cohort-{INTERVAL}.mt'
    apply_annotations(
        vcf_path=vcf_path,
        vep_ht_path=vep_ht_path,
        out_mt_path=str(out_mt_path),
        checkpoints_bucket=f'{tmp_bucket}/checkpoints',
        clinvar_ht_path=reference_path('seqr/clinvar_ht'),
    )
    # Testing
    mt = hl.read_matrix_table(str(out_mt_path))
    # mt.rows().show()
    assert mt.topmed.AC.collect() == [20555, 359, 20187]
    assert set(mt.geneIds.collect()[0]) == {'ENSG00000089063'}
