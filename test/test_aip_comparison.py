"""
test class for AIP comparisons
"""


import logging

import pytest
from peddy import Ped

import hail as hl

from comparison.comparison import (
    check_gene_is_green,
    check_variant_was_normalised,
    check_in_vcf,
    common_format_aip,
    common_format_seqr,
    find_affected_samples,
    find_missing,
    find_variant_in_mt,
    run_ac_check,
    run_quality_flag_check,
    CommonFormatResult,
    Confidence,
)


def test_parse_aip(output_json):
    """
    tests that the AIP output JSON is parsed correctly
    :param output_json:
    :return:
    """

    parsed_result = common_format_aip(output_json)
    assert list(parsed_result.keys()) == ['SAMPLE_1']

    parsed_variants = parsed_result['SAMPLE_1']
    assert isinstance(parsed_variants, list)
    assert len(parsed_variants) == 1

    parsed_var = parsed_variants[0]
    assert isinstance(parsed_var, CommonFormatResult)

    test_variant = CommonFormatResult('7', 93105286, 'T', 'A', [Confidence.EXPECTED])
    assert test_variant == parsed_var


def test_common_format_result():
    """
    test on the generic object
    :return:
    """
    common = CommonFormatResult('7', 93105286, 'T', 'A', [Confidence.EXPECTED])

    # check that the normalisation works ok
    assert (
        common
        == CommonFormatResult('chr7', 93105286, 'T', 'A', [Confidence.EXPECTED])
        == CommonFormatResult('chR7', 93105286, 'T', 'A', [Confidence.EXPECTED])
    )

    assert common.normalise_chrom('CHRchr7') == 'CHR7'


def test_affected_finder(trio_ped):
    """
    tests function to find probands from a Ped file
    this trio contains PROBAND, MOTHER, and FATHER.
    Only the proband is affected, and the mother and father both have a listed child
    :param trio_ped:
    :return:
    """
    ped_parsed = Ped(trio_ped)
    samples = find_affected_samples(ped_parsed)
    assert samples == ['PROBAND1', 'PROBAND2']


def test_affected_finder_with_sibling(quad_ped):
    """
    tests function to find probands from a Ped file
    contains the same trio as above ^^ plus an unaffected sibling
    :param quad_ped:
    """
    samples = find_affected_samples(quad_ped)
    assert samples == ['PROBAND']


def test_seqr_parser(seqr_csv_output):
    """
    First variant is 17:10697288 G>A & tagged as Possible
    Second variant does not have an AIP training tag, so should be ignored
    :param seqr_csv_output:
    :return:
    """

    seqr_results = common_format_seqr(seqr_csv_output, affected=['PROBAND'])

    # only results for one sample
    assert len(list(seqr_results.keys())) == 1

    sample_variants = seqr_results.get('PROBAND')
    assert (
        len(sample_variants) == 1
    ), f'Should only be one retained variant: {sample_variants}'

    assert sample_variants == [
        CommonFormatResult('17', 10697288, 'G', 'A', [Confidence.POSSIBLE])
    ]


def test_find_missing_matched(caplog):
    """
    trial the comparison process, check logged results
    one matching variant, and one bonus AIP result
    :param caplog:
    :return:
    """
    caplog.set_level(logging.INFO)

    fake_seqr = {'match': [CommonFormatResult('1', 1, 'A', 'C', [Confidence.POSSIBLE])]}

    fake_aip = {
        'match': [
            CommonFormatResult('1', 1, 'A', 'C', [Confidence.POSSIBLE]),
            CommonFormatResult('2', 2, 'A', 'G', [Confidence.POSSIBLE]),
        ]
    }

    discrep = find_missing(fake_aip, fake_seqr)
    assert len(discrep) == 0
    log_records = [rec.message for rec in caplog.records]
    assert 'Sample match - 1 matched variant(s)' in log_records
    assert 'Sample match - 0 missing variant(s)' in log_records


def test_find_missing_fails(caplog):
    """
    trial the comparison process, check logged results
    one mis-matching variant
    :param caplog:
    :return:
    """
    caplog.set_level(logging.INFO)

    fake_seqr = {'match': [CommonFormatResult('1', 1, 'A', 'C', [Confidence.POSSIBLE])]}

    fake_aip = {
        'match': [
            CommonFormatResult('2', 2, 'A', 'G', [Confidence.POSSIBLE]),
        ]
    }

    discrep = find_missing(fake_aip, fake_seqr)
    assert len(discrep) == 1
    assert len(discrep['match']) == 1
    log_records = [rec.message for rec in caplog.records]
    assert 'Sample match - 0 matched variant(s)' in log_records
    assert 'Sample match - 1 missing variant(s)' in log_records


def test_find_missing_different_sample(caplog):
    """
    trial the comparison process, check logged results
    no common samples
    :param caplog:
    :return:
    """
    caplog.set_level(logging.INFO)

    fake_seqr = {'match': [CommonFormatResult('1', 1, 'A', 'C', [Confidence.POSSIBLE])]}

    fake_aip = {
        'mismatch': [
            CommonFormatResult('2', 2, 'A', 'G', [Confidence.POSSIBLE]),
        ]
    }

    discrep = find_missing(fake_aip, fake_seqr)
    assert len(discrep) == 1
    assert len(discrep['match']) == 1
    log_records = [rec.message for rec in caplog.records]
    assert 'Samples completely missing from AIP results: match' in log_records
    assert 'Sample match: 1 missing variant(s)' in log_records


def test_missing_in_vcf(single_variant_vcf_path):
    """
    test method which scans VCF for a given variant
    :return:
    """
    common_var = CommonFormatResult('1', 1, 'GC', 'G', [Confidence.EXPECTED])
    variant_object = {'SAMPLE': [common_var]}
    in_vcf, not_in_vcf = check_in_vcf(single_variant_vcf_path, variant_object)
    assert len(in_vcf) == 1
    assert in_vcf['SAMPLE'] == [common_var]
    assert len(not_in_vcf) == 0


def test_missing_in_vcf_confidence_irrelevant(single_variant_vcf_path):
    """
    test method which scans VCF for a given variant
    :return:
    """
    common_var = CommonFormatResult('1', 1, 'GC', 'G', [Confidence.UNLIKELY])
    variant_object = {'SAMPLE': [common_var]}
    in_vcf, not_in_vcf = check_in_vcf(single_variant_vcf_path, variant_object)
    assert len(in_vcf) == 1
    assert in_vcf['SAMPLE'] == [common_var]
    assert len(not_in_vcf) == 0


def test_missing_in_vcf_allele_mismatch(single_variant_vcf_path):
    """
    test method which scans VCF for a given variant
    alleles should not match
    :return:
    """
    common_var = CommonFormatResult('1', 1, 'GC', 'A', [Confidence.EXPECTED])
    variant_object = {'SAMPLE': [common_var]}
    in_vcf, not_in_vcf = check_in_vcf(single_variant_vcf_path, variant_object)
    assert len(not_in_vcf) == 1
    assert not_in_vcf['SAMPLE'] == [common_var]
    assert len(in_vcf) == 0


def test_missing_in_vcf_variant_pos_mismatch(single_variant_vcf_path):
    """
    test method which scans VCF for a given variant
    alleles should not match
    :return:
    """
    common_var = CommonFormatResult('11', 11, 'GC', 'G', [Confidence.EXPECTED])
    variant_object = {'SAMPLE': [common_var]}
    in_vcf, not_in_vcf = check_in_vcf(single_variant_vcf_path, variant_object)
    assert len(not_in_vcf) == 1
    assert not_in_vcf['SAMPLE'] == [common_var]
    assert len(in_vcf) == 0


def test_variant_in_mt(hail_matrix):
    """
    check that a variant can be retrieved from a MT
    :return:
    """
    query_variant = CommonFormatResult('1', 1, 'GC', 'G', [])
    result = find_variant_in_mt(hail_matrix, query_variant)
    assert result.count_rows() == 1


def test_variant_in_mt_allele_mismatch(hail_matrix):
    """
    check that a variant can be retrieved from a MT
    :return:
    """
    query_variant = CommonFormatResult('1', 1, 'A', 'T', [])
    result = find_variant_in_mt(hail_matrix, query_variant)
    assert result.count_rows() == 0


def test_variant_in_mt_wrong_chrom(hail_matrix):
    """
    check that a variant can be retrieved from a MT
    """
    query_variant = CommonFormatResult('11', 1, 'GC', 'G', [])
    result = find_variant_in_mt(hail_matrix, query_variant)
    assert result.count_rows() == 0


def test_variant_in_mt_wrong_pos(hail_matrix):
    """
    check that a variant can be retrieved from a MT
    """
    query_variant = CommonFormatResult('1', 100, 'GC', 'G', [])
    result = find_variant_in_mt(hail_matrix, query_variant)
    assert result.count_rows() == 0


# - - - - - - - - - - - - - - - - - - - #
# tests relating to MT category/quality #
# - - - - - - - - - - - - - - - - - - - #


@pytest.mark.parametrize('gene,rows', (('green', 1), ('other', 0)))
def test_check_gene_is_green(gene, rows, hail_matrix):
    """
    test that the geneIds filter works
    """
    green_genes = hl.set(['green'])
    gene_mt = hail_matrix.annotate_rows(geneIds=gene)
    result = check_gene_is_green(gene_mt, green_genes)
    assert result.count_rows() == rows


@pytest.mark.parametrize(
    'ac,an,clinvar,results',
    [
        (100, 100, 0, ['QC: AC too high in joint call']),
        (100, 100, 1, []),
        (1, 100, 0, []),
    ],
)
def test_ac_threshold(
    ac: int, an: int, clinvar: int, results: list[str], hail_matrix: hl.MatrixTable
):
    """
    required fields: alleles, AC, AN
    """
    anno_mt = hail_matrix.annotate_rows(AC=ac, AN=an)
    anno_mt = anno_mt.annotate_rows(info=anno_mt.info.annotate(clinvar_aip=clinvar))
    assert run_ac_check(anno_mt) == results


@pytest.mark.parametrize(
    'filters,results',
    [({'FAILURE'}, ['QC: Variant has assigned quality flags']), (None, [])],
)
def test_run_quality_flag_check(
    filters: list[str],
    results: list[str],
    hail_matrix: hl.MatrixTable,
):
    """
    ronseal
    """
    if filters is None:
        filters = hl.empty_set(t=hl.tstr)
    anno_mt = hail_matrix.annotate_rows(filters=filters)
    assert run_quality_flag_check(anno_mt) == results


@pytest.mark.parametrize(
    'alleles,results',
    [
        (['A', 'C'], []),
        (['A', 'C', 'G'], ['QC: Variant not well normalised']),
        (['A', '*'], ['QC: Variant not well normalised']),
    ],
)
def test_variant_is_normalised(
    alleles: list[str],
    results: list[str],
    hail_matrix: hl.MatrixTable,
):
    """
    ronseal
    """
    anno_mt = hail_matrix.key_rows_by(hail_matrix.locus).annotate_rows(alleles=alleles)
    assert check_variant_was_normalised(anno_mt) == results
