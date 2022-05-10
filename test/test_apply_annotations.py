"""
tests for the annotation module
"""


import string
import time
from random import choices
import hail as hl
from reanalysis.annotation import apply_annotations


def timestamp():
    """
    Generate a timestamp plus a short random string for guaranteed uniqueness.
    """
    rand_bit = ''.join(choices(string.ascii_uppercase + string.digits, k=3))
    return time.strftime('%Y-%m%d-%H%M') + rand_bit


TMP_BUCKET = (
    f'gs://cpg-fewgenomes-test/unittest/tmp/test_apply_annotations/{timestamp()}'
)
INTERVAL = 'chr20-5111495-5111607'


def test_apply_annotations(monkeypatch):
    """
    Test annotation.apply_annotations: convert VCF to a MatrixTable,
    and add VEP and other annotations.
    """
    monkeypatch.setenv('CPG_REFERENCE_PREFIX', 'gs://cpg-reference')
    vcf_path = (
        f'gs://cpg-fewgenomes-test/unittest/inputs/chr20/'
        f'joint-called-{INTERVAL}.vcf.gz'
    )
    vep_ht_path = f'gs://cpg-fewgenomes-test/unittest/inputs/chr20/vep/{INTERVAL}.ht'
    out_mt_path = f'{TMP_BUCKET}/cohort-{INTERVAL}.mt'
    apply_annotations(
        vcf_path=vcf_path,
        vep_ht_path=vep_ht_path,
        out_mt_path=str(out_mt_path),
        checkpoints_bucket=f'{TMP_BUCKET}/checkpoints',
    )
    # Testing
    mt = hl.read_matrix_table(str(out_mt_path))
    # mt.rows().show()
    assert mt.topmed.AC.collect() == [20555, 359, 20187]
    assert set(mt.geneIds.collect()[0]) == {'ENSG00000089063'}
