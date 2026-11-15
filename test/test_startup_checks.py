"""
Tests for the VCF-contract startup checks (startup_checks.py).:
the annotated shard's INFO/BCSQ contract
pedigree-input data overlap
"""

import pytest
from mendelbrot.pedigree_parser import PedigreeParser

from talos import startup_checks
from talos.startup_checks import check_vcf, validate_pedigree

BCSQ_HEADER_LINE = (
    '##INFO=<ID=BCSQ,Number=.,Type=String,Description="Local consequence annotation from BCFtools/csq, '
    'Format: Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change">'
)

INFO_HEADER_LINES = {
    'AC': '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">',
    'AN': '##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">',
    'gnomad_AC': '##INFO=<ID=gnomad_AC,Number=1,Type=Integer,Description="gnomAD AC">',
    'gnomad_AF': '##INFO=<ID=gnomad_AF,Number=1,Type=Float,Description="gnomAD AF">',
    'gnomad_HomAlt': '##INFO=<ID=gnomad_HomAlt,Number=1,Type=Integer,Description="gnomAD HomAlt">',
}

BCSQ_ENTRY = 'missense|GENE1|ENST01|protein_coding|+|123P|456A>G'
BASE_INFO = 'AC=1;AN=6;gnomad_AC=1;gnomad_AF=0.0001;gnomad_AC=0;gnomad_HomAlt=0'

CLINVAR_INFO_LINES = (
    '##INFO=<ID=allele_id,Number=1,Type=Integer,Description="ClinVar Allele ID">\n'
    '##INFO=<ID=gold_stars,Number=1,Type=Integer,Description="Stars">\n'
    '##INFO=<ID=clinical_significance,Number=1,Type=String,Description="Significance">\n'
)


@pytest.fixture(autouse=True)
def _clear_errors() -> None:
    """The check functions accumulate into module-level lists, so reset them before every test."""
    startup_checks.LOG_ERRORS.clear()
    startup_checks.CONFIG_ERRORS.clear()


def write_annotated_vcf(
    tmp_path,
    *,
    rows: list[str] | None = None,
    drop_info: tuple[str, ...] = (),
    with_bcsq_header: bool = True,
    samples: tuple[str, ...] = ('proband', 'father_1', 'mother_1'),
):
    """Build a minimal stand-in for one annotated shard."""
    header_lines = [
        '##fileformat=VCFv4.2',
        '##contig=<ID=chr1,length=248956422>',
        *[line for field, line in INFO_HEADER_LINES.items() if field not in drop_info],
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    ]
    if with_bcsq_header:
        header_lines.append(BCSQ_HEADER_LINE)

    columns = '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', *samples])
    gts = '\t'.join(['0/1'] * len(samples))

    if rows is None:
        rows = [f'chr1\t12345\t.\tA\tG\t60\tPASS\tBCSQ={BCSQ_ENTRY};{BASE_INFO}\tGT\t{gts}']

    vcf_path = tmp_path / 'annotated.vcf'
    vcf_path.write_text('\n'.join([*header_lines, columns, *rows]) + '\n')
    return str(vcf_path)


def write_pedigree(tmp_path, *, sample_ids=('proband', 'father_1', 'mother_1'), affected='proband'):
    """Trio pedigree, with control over which member (if any) is affected."""
    child, father, mother = sample_ids
    rows = [
        f'family_1\t{child}\t{father}\t{mother}\t1\t{2 if affected == child else 1}',
        f'family_1\t{father}\t0\t0\t1\t{2 if affected == father else 1}',
        f'family_1\t{mother}\t0\t0\t2\t{2 if affected == mother else 1}',
    ]
    ped_path = tmp_path / 'pedigree.ped'
    ped_path.write_text('\n'.join(rows) + '\n')
    return PedigreeParser(str(ped_path))


def test_complete_vcf_passes(tmp_path):
    check_vcf(write_annotated_vcf(tmp_path), write_pedigree(tmp_path))
    assert not startup_checks.LOG_ERRORS


@pytest.mark.parametrize('missing_field', ['AC', 'AN', 'gnomad_AF', 'gnomad_HomAlt'])
def test_missing_info_field_is_reported(tmp_path, missing_field):
    check_vcf(write_annotated_vcf(tmp_path, drop_info=(missing_field,)), write_pedigree(tmp_path))
    assert any(f'INFO/{missing_field} is missing' in error for error in startup_checks.LOG_ERRORS)


def test_missing_bcsq_header_reported_without_crashing(tmp_path):
    check_vcf(write_annotated_vcf(tmp_path, with_bcsq_header=False), write_pedigree(tmp_path))
    assert any('INFO/BCSQ is missing' in error for error in startup_checks.LOG_ERRORS)


def test_bcsq_longer_than_its_header_is_reported(tmp_path):
    # seven fields are declared, this row carries nine
    overlong = f'{BCSQ_ENTRY}|extra|fields'
    row = f'chr1\t12345\t.\tA\tG\t60\tPASS\tBCSQ={overlong};{BASE_INFO}\tGT\t0/1\t0/1\t0/1'
    check_vcf(write_annotated_vcf(tmp_path, rows=[row]), write_pedigree(tmp_path))
    assert any('more fields than the 7 its header declares' in error for error in startup_checks.LOG_ERRORS)


def test_truncated_bcsq_is_accepted(tmp_path):
    # BCFtools drops trailing empty fields, so a short entry is normal
    row = f'chr1\t12345\t.\tA\tG\t60\tPASS\tBCSQ=missense|GENE1|ENST01;{BASE_INFO}\tGT\t0/1\t0/1\t0/1'
    check_vcf(write_annotated_vcf(tmp_path, rows=[row]), write_pedigree(tmp_path))
    assert not startup_checks.LOG_ERRORS


def test_empty_vcf_is_reported(tmp_path, caplog):
    check_vcf(write_annotated_vcf(tmp_path, rows=[]), write_pedigree(tmp_path))
    assert 'contains no variants' in caplog.text


def test_rows_without_bcsq_are_reported(tmp_path, caplog):
    row = f'chr1\t12345\t.\tA\tG\t60\tPASS\t{BASE_INFO}\tGT\t0/1\t0/1\t0/1'
    check_vcf(write_annotated_vcf(tmp_path, rows=[row]), write_pedigree(tmp_path))
    assert 'carry a BCSQ annotation' in caplog.text


def test_pedigree_with_no_shared_samples_is_reported(tmp_path):
    pedigree = write_pedigree(tmp_path, sample_ids=('other_1', 'other_2', 'other_3'), affected='other_1')
    check_vcf(write_annotated_vcf(tmp_path), pedigree)
    assert any('No pedigree samples are present' in error for error in startup_checks.LOG_ERRORS)


def test_pedigree_with_no_affected_sample_in_vcf_is_reported(tmp_path):
    # the affected member is absent from the VCF, the unaffected parents are present
    pedigree = write_pedigree(tmp_path, sample_ids=('absent_child', 'father_1', 'mother_1'), affected='absent_child')
    check_vcf(write_annotated_vcf(tmp_path), pedigree)
    assert any('No affected pedigree samples are present' in error for error in startup_checks.LOG_ERRORS)


def test_unreadable_pedigree_is_reported(tmp_path):
    assert validate_pedigree(str(tmp_path / 'no_such_file.ped')) is None
    assert any('Error parsing pedigree file' in error for error in startup_checks.LOG_ERRORS)


def test_pedigree_without_affected_members_is_reported(tmp_path):
    write_pedigree(tmp_path, affected='nobody')
    validate_pedigree(str(tmp_path / 'pedigree.ped'))
    assert any('does not contain affected members' in error for error in startup_checks.LOG_ERRORS)
