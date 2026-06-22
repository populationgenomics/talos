"""
code for lifting over models from 1.1.0 to 1.2.0
"""


def resultdata(data_dict: dict) -> dict:
    """
    Lift over ResultData from 1.1.0 to 1.2.0
    Requires the adjustment of ResultData.results.ParticipantResults.metadata.ParticipantMeta.panel_*
    """
    # check we're upgrading the right version
    # this could be from any prior version
    if not data_dict['version'] < '1.2.0':
        raise AssertionError(f'This method cannot upgrade from {data_dict["version"]}')

    for _sample, content in data_dict['results'].items():
        for variant in content['variants']:
            # this attribute isn't exported by default, so allow it to be missing
            if 'sample_support' in variant:
                _ = variant['var_data'].pop('sample_support')

            # these fields did not exist before 1.2.0; default them to empty rather than reading live config.
            # a historical document carries no record of the support/ignore categories used at the time, and a
            # migration must be a pure function of the old document (not of whatever config happens to be loaded now)
            variant['var_data']['support_categories'] = []
            variant['var_data']['ignored_categories'] = []

    data_dict['version'] = '1.2.0'
    return data_dict
