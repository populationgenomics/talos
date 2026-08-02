"""
Integration tests for the SVAFotate step of the SV_ANNOTATION workflow.

These invoke the svafotate CLI, which lives only in the SVAFotate container (its dependencies are
irreconcilable with Talos's - see docs/SvAnnotationPlan.md). They are skipped when svafotate is not on PATH,
so they no-op in unit CI and run inside the SVAFotate image / an integration environment.

They guard SVAFotate behaviours the plan flags as silent-failure traps, which the pure-Python unit tests in
test_rename_sv_af_fields.py cannot reach:

- the shipped Ensembl-style reference BED annotates a chr-prefixed query VCF with non-zero AFs, while a
  chr-prefixed (contig-normalised) BED silently zeroes every frequency. This is the single most important
  guard - "normalising" the BED would zero an entire callset with no error raised.
- the reciprocal-overlap threshold (-f) stops a small variant inheriting a large record's frequency
- end to end with the rename, a sub-threshold candidate overlap yields a null SVID, not the spurious
  Best_gnomAD_ID that SVAFotate reports regardless of the threshold

The reference is test/input/svafotate_slice.bed - a two-record slice of the real
SVAFotate_SV_popAFs.GRCh38.v4.1 BED (one DUP, one DEL), kept Ensembl-style (contig '1') exactly as shipped.
"""

import gzip
import shutil
import subprocess
from pathlib import Path

import pytest
from cyvcf2 import VCFReader

from talos.annotation_scripts import rename_sv_af_fields

SVAFOTATE = shutil.which('svafotate')
pytestmark = pytest.mark.skipif(SVAFOTATE is None, reason='svafotate CLI not on PATH')

BED_SLICE = Path(__file__).parent / 'input' / 'svafotate_slice.bed'
OVERLAP_FRACTION = '0.5'

# the DUP record in the slice: gnomAD-SV_v3_DUP_chr1_a5cef17b, 1:63999-74000, AF ~0.1188
DUP_AF = 0.118818
DUP_SVID = 'gnomAD-SV_v3_DUP_chr1_a5cef17b'
# the DEL record in the slice: gnomAD-SV_v3_DEL_chr1_143f6519, 1:509967-690000 (180kb), AF ~0.2238

QUERY_VCF = """\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=DEL,Description="Deletion">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample
chr1\t64000\tdup_match\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=74000;SVLEN=10000\tGT\t0/1
chr1\t600000\tdel_small\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=600200;SVLEN=200\tGT\t0/1
"""


def _write_query_vcf(tmp_path) -> Path:
    """the chr-prefixed query: a DUP with ~100% overlap of the slice's DUP, and a 200bp DEL deep inside the
    180kb DEL (far below -f 0.5)"""
    path = tmp_path / 'query.vcf'
    path.write_text(QUERY_VCF)
    return path


def _write_bed(tmp_path, add_chr: bool) -> Path:
    """
    write the slice as a bgzipped BED, optionally rewriting the data rows to UCSC-style contigs. The header
    line is left untouched; only the leading contig of each data row is prefixed
    """
    out_lines = []
    for line in BED_SLICE.read_text().splitlines():
        out_lines.append(f'chr{line}' if add_chr and not line.startswith('#') else line)
    dest = tmp_path / ('bed_chr.bed.gz' if add_chr else 'bed_ensembl.bed.gz')
    with gzip.open(dest, 'wt') as handle:
        handle.write('\n'.join(out_lines) + '\n')
    return dest


def _run_svafotate(tmp_path, bed: Path) -> Path:
    """run svafotate with the flags the AnnotateSvWithSvafotate module uses, return the output VCF path"""
    out = tmp_path / 'popaf.vcf'
    subprocess.run(  # noqa: S603
        [
            SVAFOTATE,
            'annotate',
            '-v',
            str(_write_query_vcf(tmp_path)),
            '-o',
            str(out),
            '-b',
            str(bed),
            '-s',
            'gnomAD',
            '-f',
            OVERLAP_FRACTION,
            '-a',
            'best',
            '--ins',
            '--cpu',
            '1',
        ],
        check=True,
    )
    return out


def _by_id(vcf_path: Path) -> dict:
    """parse a VCF into an ID -> variant lookup"""
    return {variant.ID: variant for variant in VCFReader(str(vcf_path))}


def test_ensembl_bed_annotates_chr_query(tmp_path):
    """the shipped Ensembl-style BED annotates a chr-prefixed query with the real gnomAD frequency"""
    variants = _by_id(_run_svafotate(tmp_path, _write_bed(tmp_path, add_chr=False)))
    assert variants['dup_match'].INFO.get('Max_AF') == pytest.approx(DUP_AF, rel=1e-3)


def test_chr_prefixed_bed_silently_zeroes(tmp_path):
    """
    the silent-failure guard: a contig-normalised (chr-prefixed) BED matches nothing, so every frequency
    reads 0 with no error. This is what would happen if someone "fixed" the BED's contigs
    """
    variants = _by_id(_run_svafotate(tmp_path, _write_bed(tmp_path, add_chr=True)))
    assert variants['dup_match'].INFO.get('Max_AF') == 0


def test_overlap_threshold_blocks_small_variant(tmp_path):
    """a 200bp DEL deep inside a 180kb DEL is below -f 0.5, so it must not inherit the large record's AF"""
    variants = _by_id(_run_svafotate(tmp_path, _write_bed(tmp_path, add_chr=False)))
    assert variants['del_small'].INFO.get('Max_AF') == 0


def test_subthreshold_svid_withheld_after_rename(tmp_path):
    """
    end to end with the rename: the matching DUP keeps its SVID, but the sub-threshold DEL - for which
    SVAFotate still reports a spurious Best_gnomAD_ID - is emitted with a null SVID
    """
    popaf = _run_svafotate(tmp_path, _write_bed(tmp_path, add_chr=False))
    renamed = tmp_path / 'renamed.vcf'
    rename_sv_af_fields.main(vcf_in=str(popaf), vcf_out=str(renamed))
    variants = _by_id(renamed)

    assert variants['dup_match'].INFO.get(rename_sv_af_fields.TALOS_AF) == pytest.approx(DUP_AF, rel=1e-3)
    assert variants['dup_match'].INFO.get(rename_sv_af_fields.TALOS_SVID) == DUP_SVID
    assert variants['del_small'].INFO.get(rename_sv_af_fields.TALOS_SVID) is None
