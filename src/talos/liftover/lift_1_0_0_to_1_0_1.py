"""
code for lifting over models from 1.0.0 to 1.0.1
"""


def resultdata(data_dict: dict) -> dict:
    """
    Lift over ResultData from 1.0.0 to 1.0.1
    Requires the adjustment of "first_seen" to "first_tagged"
    """
    # check we're upgrading the right version
    # this could be from any prior version
    assert data_dict['version'] < '1.0.1'
    for results in data_dict['results'].values():
        for rv in results['variants']:
            rv['date_of_phenotype_match'] = None
            first = rv.pop('first_seen')
            rv['first_tagged'] = first
            rv['evidence_last_updated'] = first

    return data_dict


def historicvariants(data_dict: dict) -> dict:
    """
    Lift over HistoricVariants from 1.0.0 to 1.0.1
    Requires the adjustment of "first_seen" to "first_tagged"
    in results.sample.variant.first_seen
    """
    assert data_dict['version'] < '1.0.1'
    for results in data_dict['results'].values():
        for details in results.values():
            earliest = min(details['categories'].values())
            details['first_tagged'] = earliest
            # this doesn't matter, pydantic can handle it
            details['support_vars'] = set(details.get('support_vars', []))
    return data_dict
