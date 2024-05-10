"""
code for lifting over models from None/unspecified to 1.0.0
"""


def lift_pmp(data_dict: dict) -> dict:
    """
    Lift over PhenotypeMatchedPanels from None to 1.0.0
    requires the migration of HPO terms from strings to dicts
    prev: ['HPO:0000001 - Definition', ]
    required: [{'id': 'HPO:0000001', 'label': 'Definition'}, ]
    """
    # confirm that we're upgrading the right version
    assert data_dict.get('version') is None
    for _sample, data in data_dict['samples'].items():
        assert isinstance(data['hpo_terms'], list)
        new_data: list[dict] = []
        for term in data['hpo_terms']:
            try:
                # expect two values
                hpo_id, label = term.split(' - ')
                new_data.append({'id': hpo_id, 'label': label})

            except IndexError:
                assert isinstance(term, str)
                hpo_id = term
                new_data.append({'id': hpo_id, 'label': ''})

            except AttributeError:
                # wasn't a string?
                raise ValueError(f'Unexpected HPO term format: {term}')

        data['hpo_terms'] = new_data
    return data_dict


def lift_resultdata(data_dict: dict) -> dict:
    """
    Lift over ResultData from None to 1.0.0
    requires the migration of HPO terms from strings to dicts
    prev: ['HPO:0000001 - Definition', ]
    required: [{'id': 'HPO:0000001', 'label': 'Definition'}, ]
    """
    # confirm that we're upgrading the right version
    assert data_dict.get('version') is None
    for _sample, data in data_dict['results'].items():
        assert isinstance(data['metadata']['phenotypes'], list)
        new_data: list[dict] = []
        for term in data['metadata']['phenotypes']:
            try:
                # expect two values
                hpo_id, label = term.split(' - ')
                new_data.append({'id': hpo_id, 'label': label})

            except ValueError:
                assert isinstance(term, str)
                hpo_id = term
                new_data.append({'id': hpo_id, 'label': ''})

            except AttributeError:
                # wasn't a string?
                raise ValueError(f'Unexpected HPO term format: {term}')

        data['metadata']['phenotypes'] = new_data
    return data_dict
