"""
test class for AIP comparisons
"""


import logging
from peddy import Ped

from reanalysis.comparison import (
    common_format_from_results,
    common_format_from_seqr,
    find_affected_samples,
    find_missing,
    CommonFormatResult,
    Confidence,
)


def test_parse_aip(output_json):
    """
    tests that the AIP output JSON is parsed correctly
    :param output_json:
    :return:
    """

    parsed_result = common_format_from_results(output_json)
    assert list(parsed_result.keys()) == ['SAMPLE_1']

    parsed_variants = parsed_result['SAMPLE_1']
    assert isinstance(parsed_variants, list)
    assert len(parsed_variants) == 1

    parsed_var = parsed_variants[0]
    assert isinstance(parsed_var, CommonFormatResult)

    test_variant = CommonFormatResult('7', 93105286, 'T', 'A', [Confidence.EXPECTED])
    assert test_variant == parsed_var


def test_affected_finder(trio_ped):
    """
    tests function to find probands from a Ped file
    this trio contains PROBAND, MOTHER, and FATHER.
    Only the proband is affected, and the mother and father both have a listed child
    not the toughest test
    :param trio_ped:
    :return:
    """
    # digest that Ped
    ped_parsed = Ped(trio_ped)
    samples = find_affected_samples(ped_parsed)
    assert samples == ['PROBAND']


def test_affected_finder_with_sibling(quad_ped):
    """
    tests function to find probands from a Ped file
    contains the same trio as above ^^ plus an unaffected sibling
    :param quad_ped:
    :return:
    """
    # digest that Ped
    ped_parsed = Ped(quad_ped)
    samples = find_affected_samples(ped_parsed)
    assert samples == ['PROBAND']


def test_seqr_parser(seqr_csv_output):
    """
    First variant is 17:10697288 G>A & tagged as Possible
    Second variant does not have an AIP training tag, so should be ignored
    :param seqr_csv_output:
    :return:
    """

    seqr_results = common_format_from_seqr(seqr_csv_output, affected=['PROBAND'])

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
