"""
Tests for talos.scripts.prepare_sv_for_svannotate - the pass that prepares an SV VCF for GATK SVAnnotate.

Two behaviours are covered:
- records with any non-symbolic ALT are removed (SVAnnotate only works with symbolic alleles), and each is logged
- CPX_INTERVALS is coordinate-sorted on the records that remain
"""

import pytest
from cyvcf2 import VCFReader
from loguru import logger

from talos.scripts.prepare_sv_for_svannotate import is_symbolic, main, sorted_intervals

VCF_HEADER = """\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##ALT=<ID=DEL,Description="Deletion">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">
##INFO=<ID=CPX_INTERVALS,Number=.,Type=String,Description="Complex SV intervals">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


@pytest.mark.parametrize(
    'allele',
    ['<DEL>', '<INS:ME:ALU>', '<CPX>', 'N[chr2:123[', ']chr2:123]N', 'N.', '.N'],
)
def test_is_symbolic(allele):
    """bracketed, breakpoint and single-breakend alleles are all symbolic to htsjdk"""
    assert is_symbolic(allele)


@pytest.mark.parametrize('allele', ['A', 'ACGTACGT', '*', '.', ''])
def test_is_not_symbolic(allele):
    """sequence-resolved alleles, single bases and placeholders are not symbolic"""
    assert not is_symbolic(allele)


def test_sorted_intervals_reorders_delinvdel():
    """GATK-SV's DEL,DEL,INV ordering is put into coordinate order"""
    raw = 'DEL_chr1:1416597-16725244,DEL_chr1:64416883-234776440,INV_chr1:16725244-64416883'
    assert sorted_intervals(raw) == 'DEL_chr1:1416597-16725244,INV_chr1:16725244-64416883,DEL_chr1:64416883-234776440'


def _run(tmp_path, rows: list[str]) -> tuple[list, list[str]]:
    """run the pass over a fixture VCF, return the output variants and the captured log messages"""
    vcf_in = tmp_path / 'in.vcf'
    vcf_in.write_text(VCF_HEADER + '\n'.join(rows) + '\n')
    vcf_out = str(tmp_path / 'out.vcf')

    messages: list[str] = []
    handler = logger.add(messages.append, format='{message}')
    try:
        main(vcf_in=str(vcf_in), vcf_out=vcf_out)
    finally:
        logger.remove(handler)
    return list(VCFReader(vcf_out)), messages


def test_non_symbolic_records_removed_and_logged(tmp_path):
    """only fully-symbolic records survive, and every removed record is logged with its site columns"""
    variants, messages = _run(
        tmp_path,
        [
            'chr1\t1000\tkeep_del\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=5000',
            'chr1\t2000\tdrop_seq\tA\tACGTACGTACGT\t.\tPASS\tSVTYPE=INS;END=2001',
            'chr1\t3000\tkeep_bnd\tN\tN[chr2:123[\t.\tPASS\tSVTYPE=BND;END=3000',
            'chr1\t4000\tdrop_mixed\tN\t<DEL>,ACGT\t.\tPASS\tSVTYPE=DEL;END=6000',
            'chr1\t5000\tdrop_no_alt\tN\t.\t.\tPASS\tSVTYPE=DEL;END=7000',
        ],
    )

    assert [variant.ID for variant in variants] == ['keep_del', 'keep_bnd']

    removed = [message for message in messages if message.startswith('Removing non-symbolic record')]
    assert len(removed) == 3
    assert 'chr1\t2000\tdrop_seq\tA\tACGTACGTACGT' in removed[0]
    assert 'chr1\t4000\tdrop_mixed\tN\t<DEL>,ACGT' in removed[1]
    assert 'chr1\t5000\tdrop_no_alt\tN\t.' in removed[2]
    assert any('Removed 3 of 5 records' in message for message in messages)


def test_cpx_intervals_sorted_on_kept_records(tmp_path):
    """the CPX_INTERVALS sort still applies to records that pass the symbolic filter"""
    (variant,) = _run(
        tmp_path,
        [
            'chr1\t1416597\tcpx\tN\t<CPX>\t.\tPASS\tSVTYPE=CPX;END=234776440;'
            'CPX_INTERVALS=DEL_chr1:1416597-16725244,DEL_chr1:64416883-234776440,INV_chr1:16725244-64416883',
        ],
    )[0]
    assert variant.INFO.get('CPX_INTERVALS') == (
        'DEL_chr1:1416597-16725244,INV_chr1:16725244-64416883,DEL_chr1:64416883-234776440'
    )


def test_svlen_backfilled_where_derivable(tmp_path):
    """a missing SVLEN becomes END - POS; an existing one is kept; one with END <= POS is left absent"""
    variants, _messages = _run(
        tmp_path,
        [
            'chr1\t1000\tno_len\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=5000',
            'chr1\t2000\thas_len\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=5000;SVLEN=7',
            'chr1\t3000\tbnd\tN\tN[chr2:123[\t.\tPASS\tSVTYPE=BND;END=3000',
        ],
    )
    assert [variant.INFO.get('SVLEN') for variant in variants] == [4000, 7, None]


def test_svlen_header_added_when_absent(tmp_path):
    """a callset with no SVLEN header gets one, so the backfilled value is declared"""
    vcf_in = tmp_path / 'in.vcf'
    vcf_in.write_text(
        VCF_HEADER.replace('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">\n', '')
        + 'chr1\t1000\tno_len\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=5000\n'
    )
    vcf_out = str(tmp_path / 'out.vcf')
    main(vcf_in=str(vcf_in), vcf_out=vcf_out)
    (variant,) = list(VCFReader(vcf_out))
    assert variant.INFO.get('SVLEN') == 4000
