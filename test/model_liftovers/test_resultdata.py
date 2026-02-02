"""
test that bad data throws validation errors
test that model liftover works
"""

import json
from os.path import join

from talos.models import CURRENT_VERSION, HpoTerm, ResultData, lift_up_model_version


def test_rd_from_none(test_input_models_path):
    """
    test that the lift over from None to 1.0.0 works
    """
    # load the data
    with open(join(test_input_models_path, 'result_data_v_none.json'), encoding='utf-8') as handle:
        data = json.load(handle)
    assert data.get('version') is None

    # lift it
    lifted = lift_up_model_version(data, model=ResultData)

    # check the version
    assert lifted['version'] == CURRENT_VERSION

    # ingest
    parsed = ResultData.model_validate(lifted)
    assert parsed.version == CURRENT_VERSION

    # check the HPO terms
    for _sample, results in parsed.results.items():
        for term in results.metadata.phenotypes:
            assert isinstance(term, HpoTerm)
            assert isinstance(term.id, str)
            assert isinstance(term.label, str)
            if term.id == 'HP:0000001':
                assert term.label == 'translation text'
            if term.id == 'HP:0000002':
                assert term.label == ''


def test_rd_from_1_0_0(test_input_models_path):
    """
    test that the lift over from 1.0.0 works
    """
    # load the data
    with open(join(test_input_models_path, 'result_data_v_1_0_0.json'), encoding='utf-8') as handle:
        data = json.load(handle)
    assert data.get('version') == '1.0.0'
    var_first_seen = data['results']['sam2']['variants'][0]['first_seen']

    # lift it
    lifted = lift_up_model_version(data, model=ResultData)

    # check the version
    assert lifted['version'] == CURRENT_VERSION

    # ingest
    parsed = ResultData.model_validate(lifted)
    assert parsed.version == CURRENT_VERSION

    # check the HPO terms
    for _sample, results in parsed.results.items():
        for var in results.variants:
            assert var.date_of_phenotype_match is None
            assert var.evidence_last_updated == var.first_tagged == var_first_seen
