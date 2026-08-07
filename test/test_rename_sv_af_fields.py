"""
Regression tests for talos.annotation_scripts.rename_sv_af_fields.

These guard the two field-naming rules from docs/SvAnnotationPlan.md that are easy to get wrong:
- Max_AF (the conservative maximum across matches), not Best_gnomAD_AF, is copied into {prefix}_sv_AF
- Best_gnomAD_ID is only carried through when gnomAD_Count > 0. SVAFotate populates it even for candidate
  overlaps that never cleared the -f threshold, so a naive rename stamps a spurious SVID onto nearly every
  variant. This is the single most important behaviour in the script.

No SVAFotate invocation is needed - the input is a hand-written VCF shaped like SVAFotate's output.
"""

import pytest
from cyvcf2 import VCFReader

from talos.annotation_scripts.rename_sv_af_fields import TALOS_AF, TALOS_SVID, main

VCF_HEADER = """\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=Max_AF,Number=1,Type=Float,Description="SVAFotate max AF">
##INFO=<ID=Best_gnomAD_AF,Number=1,Type=Float,Description="SVAFotate best-match AF">
##INFO=<ID=gnomAD_Count,Number=1,Type=Integer,Description="SVAFotate gnomAD match count">
##INFO=<ID=Best_gnomAD_ID,Number=1,Type=String,Description="SVAFotate best gnomAD ID">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def _write_vcf(tmp_path, info_lines: list[str]) -> str:
    """write a small SVAFotate-style VCF, one variant per INFO string, return its path"""
    rows = [f'chr1\t{1000 * (idx + 1)}\tv{idx}\tN\t<DEL>\t.\tPASS\t{info}' for idx, info in enumerate(info_lines)]
    path = tmp_path / 'svafotate_style.vcf'
    path.write_text(VCF_HEADER + '\n'.join(rows) + '\n')
    return str(path)


def _run(tmp_path, info_lines: list[str]) -> list:
    """run the rename over a fixture VCF, return the output variants as a list"""
    vcf_in = _write_vcf(tmp_path, info_lines)
    vcf_out = str(tmp_path / 'renamed.vcf')
    main(vcf_in=vcf_in, vcf_out=vcf_out)
    return list(VCFReader(vcf_out))


def test_match_copies_af_and_svid(tmp_path):
    """a qualifying match copies Max_AF and Best_gnomAD_ID into the Talos field names"""
    (variant,) = _run(
        tmp_path,
        ['SVTYPE=DEL;END=5000;Max_AF=0.1188;gnomAD_Count=2;Best_gnomAD_ID=gnomAD-SV_real'],
    )
    assert variant.INFO.get(TALOS_AF) == pytest.approx(0.1188)
    assert variant.INFO.get(TALOS_SVID) == 'gnomAD-SV_real'


def test_subthreshold_match_withholds_spurious_svid(tmp_path):
    """
    gnomAD_Count == 0 means no overlap cleared the threshold. The AF (0) is still written, but the spurious
    Best_gnomAD_ID must NOT be carried into the SVID field
    """
    (variant,) = _run(
        tmp_path,
        ['SVTYPE=DEL;END=20200;Max_AF=0;gnomAD_Count=0;Best_gnomAD_ID=gnomAD-SV_spurious'],
    )
    assert variant.INFO.get(TALOS_AF) == 0
    assert variant.INFO.get(TALOS_SVID) is None


def test_copies_max_af_not_best_af(tmp_path):
    """Max_AF (the conservative maximum), not Best_gnomAD_AF, is the value copied into {prefix}_sv_AF"""
    (variant,) = _run(
        tmp_path,
        ['SVTYPE=DEL;END=5000;Max_AF=0.25;Best_gnomAD_AF=0.05;gnomAD_Count=3;Best_gnomAD_ID=gnomAD-SV_real'],
    )
    assert variant.INFO.get(TALOS_AF) == pytest.approx(0.25)


def test_missing_max_af_defaults_to_zero(tmp_path):
    """a variant SVAFotate left without Max_AF is treated as frequency 0, not dropped or errored"""
    (variant,) = _run(tmp_path, ['SVTYPE=DEL;END=5000;gnomAD_Count=0'])
    assert variant.INFO.get(TALOS_AF) == 0
    assert variant.INFO.get(TALOS_SVID) is None


def test_all_zero_af_warns(tmp_path, caplog):
    """
    a callset in which nothing matched warns - the likely cause is a contig-normalised reference BED silently
    zeroing every frequency, the failure mode the plan calls out as needing an explicit guard
    """
    _run(tmp_path, ['SVTYPE=DEL;END=20200;Max_AF=0;gnomAD_Count=0'])
    assert 'No variant carried a non-zero Max_AF' in caplog.text


def test_non_zero_af_does_not_warn(tmp_path, caplog):
    """the opposite - when at least one variant matched, no false alarm is raised"""
    _run(
        tmp_path,
        ['SVTYPE=DEL;END=5000;Max_AF=0.1188;gnomAD_Count=2;Best_gnomAD_ID=gnomAD-SV_real'],
    )
    assert 'No variant carried a non-zero Max_AF' not in caplog.text
