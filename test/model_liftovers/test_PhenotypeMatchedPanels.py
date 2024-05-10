"""
test that bad data throws validation errors
test that model liftover works
"""

import json

import pytest

from reanalysis.models import CURRENT_VERSION, PanelApp, PhenoPacketHpo, PhenotypeMatchedPanels, lift_up_model_version
from reanalysis.utils import read_json_from_path


def test_lift_pmp_from_none(test_input_models_path):
    """
    test that the lift over from None to 1.0.0 works
    """
    # load the data
    file_path = test_input_models_path / 'hpo_panel_version_none.json'
    data = json.load(open(file_path))
    assert data.get('version') is None

    # lift it
    lifted = lift_up_model_version(data, model=PhenotypeMatchedPanels)  # type: ignore

    # check the version
    assert lifted['version'] == CURRENT_VERSION

    # ingest
    parsed = PhenotypeMatchedPanels.model_validate(lifted)

    # check the HPO terms
    for sample in parsed.samples.values():
        for term in sample.hpo_terms:
            assert isinstance(term, PhenoPacketHpo)
            assert isinstance(term.id, str)
            assert isinstance(term.label, str)
            if term.id == 'HP:0025386':
                assert term.label == 'Bitemporal hollowing'


def test_lift_pmp_from_none_json(test_input_models_path):
    """
    test that the lift over from None to 1.0.0 works
    """
    # load the data
    file_path = test_input_models_path / 'hpo_panel_version_none.json'
    lifted = read_json_from_path(file_path, return_model=PhenotypeMatchedPanels)  # type: ignore

    # check the version
    assert lifted.version == CURRENT_VERSION

    # check the HPO terms
    for sample in lifted.samples.values():
        for term in sample.hpo_terms:
            print(type(term))
            assert isinstance(term, PhenoPacketHpo)
            assert isinstance(term.id, str)
            assert isinstance(term.label, str)


def test_liftover_fails_potato(test_input_models_path):

    file_path = test_input_models_path / 'hpo_panel_version_potato.json'
    with pytest.raises(ValueError) as potato:
        read_json_from_path(file_path, return_model=PhenotypeMatchedPanels)  # type: ignore
    assert "Unknown PhenotypeMatchedPanels version: potato" in str(potato.value)


def test_no_liftover_passes(test_input_models_path, caplog):
    """test that the liftover method doesn't fail"""

    file_path = test_input_models_path / 'hpo_panel_version_none.json'
    data = json.load(open(file_path))
    assert data.get('version') is None

    # lift it
    _lifted = lift_up_model_version(data, model=PanelApp)  # type: ignore
    assert 'No liftover methods for PanelApp' in caplog.text
