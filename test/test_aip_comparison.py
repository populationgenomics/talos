"""
test class for AIP comparisons
"""

from peddy import Ped

from reanalysis.comparison import (
    common_format_from_results,
    common_format_from_seqr,
    find_probands,
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


def test_proband_finder(trio_ped):
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
    probands = find_probands(ped_parsed)
    assert probands == ['PROBAND']


def test_proband_finder_with_sibling(quad_ped):
    """
    tests function to find probands from a Ped file
    contains the same trio as above ^^ plus an unaffected sibling
    :param quad_ped:
    :return:
    """
    # digest that Ped
    ped_parsed = Ped(quad_ped)
    probands = find_probands(ped_parsed)
    assert probands == ['PROBAND']


def test_seqr_parser(seqr_csv_output):
    """
    First variant is 17:10697288 G>A & tagged as Possible
    Second variant does not have an AIP training tag, so should be ignored
    :param seqr_csv_output:
    :return:
    """

    seqr_results = common_format_from_seqr(seqr_csv_output, probands=['PROBAND'])

    # only results for one sample
    assert len(list(seqr_results.keys())) == 1

    sample_variants = seqr_results.get('PROBAND')
    assert (
        len(sample_variants) == 1
    ), f'Should only be one retained variant: {sample_variants}'

    assert sample_variants == [
        CommonFormatResult('17', 10697288, 'G', 'A', [Confidence.POSSIBLE])
    ]
