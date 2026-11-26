"""
tests for the report-building helpers in talos.create_talos_html
"""

import pytest

from talos.create_talos_html import gnomad_sv_browser_id


@pytest.mark.parametrize(
    ('svid', 'expected'),
    [
        ('gnomAD-SV_v3_DEL_chr1_143f6519', 'DEL_chr1_143f6519'),
        ('gnomAD-SV_v4_DUP_chrX_eb9a525a', 'DUP_chrX_eb9a525a'),
        # already in browser form
        ('DEL_chr1_143f6519', 'DEL_chr1_143f6519'),
    ],
)
def test_gnomad_sv_browser_id(svid: str, expected: str):
    """SVAFotate's gnomAD SV IDs are stripped to the ID the gnomAD browser resolves (dataset=gnomad_sv_r4)"""
    assert gnomad_sv_browser_id(svid) == expected
