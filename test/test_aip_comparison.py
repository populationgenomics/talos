"""
test class for AIP comparisons
"""


from reanalysis.comparison import (
    common_format_from_results,
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
