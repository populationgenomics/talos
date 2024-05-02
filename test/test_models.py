# test that bad data throws validation errors

# test that model liftover works

# not testing pydantic here, just the implementation
import json

import pytest

from reanalysis.models import PhenotypeMatchedPanels, ResultData, lift_up_model_version
from reanalysis.utils import read_json_from_path


def test_lift_pmp_none_to_1(hpo_panel_vnone):
    """
    test that the lift over from None to 1.0.0 works
    """
    # load the data
    data = json.load(open(hpo_panel_vnone))
    assert data.get('version') is None

    # lift it
    lifted = lift_up_model_version(data, model=PhenotypeMatchedPanels)

    # check the version
    assert lifted['version'] == '1.0.0'

    # ingest
    parsed = PhenotypeMatchedPanels.model_validate(lifted)

    # check the HPO terms
    for sample in parsed.samples.values():
        for term in sample.hpo_terms:
            assert isinstance(term, dict)
            assert 'id' in term
            assert 'label' in term
            assert isinstance(term['id'], str)
            assert isinstance(term['label'], str)
            if term['id'] == 'HP:0025386':
                assert term['label'] == 'Bitemporal hollowing'


def test_lift_pmp_none_to_1_from_json(hpo_panel_vnone):
    """
    test that the lift over from None to 1.0.0 works
    """
    # load the data
    data = json.load(open(hpo_panel_vnone))
    assert data.get('version') is None

    # lift it
    lifted = read_json_from_path(hpo_panel_vnone, return_model=PhenotypeMatchedPanels)

    # check the version
    assert lifted.version == '1.0.0'

    # check the HPO terms
    for sample in lifted.samples.values():
        for term in sample.hpo_terms:
            assert isinstance(term, dict)
            assert 'id' in term
            assert 'label' in term
            assert isinstance(term['id'], str)
            assert isinstance(term['label'], str)


def test_liftover_fails_potato(hpo_panel_vpotato):

    with pytest.raises(ValueError) as potato:
        read_json_from_path(hpo_panel_vpotato, return_model=PhenotypeMatchedPanels)
    assert "Unknown PhenotypeMatchedPanels version: potato" in str(potato.value)


def test_no_liftover_fails(hpo_panel_vnone):

    with pytest.raises(ValueError) as no_liftover:
        read_json_from_path(hpo_panel_vnone, return_model=ResultData)
    assert "No liftover methods for ResultData" in str(no_liftover.value)
