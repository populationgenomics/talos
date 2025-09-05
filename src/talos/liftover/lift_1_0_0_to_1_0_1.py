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
    if not data_dict['version'] < '1.0.1':
        raise AssertionError(f'This method cannot upgrade from {data_dict["version"]}')
    for results in data_dict['results'].values():
        for rv in results['variants']:
            rv['date_of_phenotype_match'] = None
            first = rv.pop('first_seen')
            rv['first_tagged'] = first
            rv['evidence_last_updated'] = first

    return data_dict
