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


def lift_fixture(test_input_models_path: str, name: str) -> dict:
    with open(join(test_input_models_path, name), encoding='utf-8') as handle:
        return lift_up_model_version(json.load(handle), model=ResultData)


def test_rd_sparse_1_0_0_lifts_like_full(test_input_models_path):
    """
    history without `independent`, with categories on var_data rather than the variant, and with
    var_data.sample_support, lifts all the way up - and ends with the same categories as the full fixture
    """
    full = lift_fixture(test_input_models_path, 'result_data_v_1_0_0.json')
    sparse = lift_fixture(test_input_models_path, 'result_data_v_1_0_0_sparse.json')
    assert sparse['version'] == CURRENT_VERSION

    full_var = full['results']['sam2']['variants'][0]
    sparse_var = sparse['results']['sam2']['variants'][0]
    assert sparse_var['categories'] == full_var['categories'] == {'2': '2024-05-09', '1': '2024-05-09'}
    assert 'independent' not in sparse_var
    assert 'sample_support' not in sparse_var['var_data']

    ResultData.model_validate(sparse)


def test_rd_2_5_0_sets_newly_green_date(test_input_models_path):
    """the 2.4.0 -> 2.5.0 placeholder lands on the model's field name"""
    variant = lift_fixture(test_input_models_path, 'result_data_v_1_0_0.json')['results']['sam2']['variants'][0]
    assert variant['newly_green_date'] == ''
    assert 'newly_green' not in variant


def test_rd_2_1_0_reasons_list_to_str(test_input_models_path):
    """the 2.0.0 -> 2.1.0 reasons list -> str step: one reason is kept, an empty list becomes ''"""
    one_reason = lift_fixture(test_input_models_path, 'result_data_v_1_0_0.json')
    assert one_reason['results']['sam2']['variants'][0]['reasons'] == 'Automsomal Dominant'

    # result_data_v1_0_1.json carries "reasons": []
    no_reason = lift_fixture(test_input_models_path, 'result_data_v1_0_1.json')
    assert no_reason['results']['sam2']['variants'][0]['reasons'] == ''
    ResultData.model_validate(no_reason)
