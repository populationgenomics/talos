"""
code for lifting over models from 1.1.0 to 1.2.0
"""

from talos.config import config_retrieve


def resultdata(data_dict: dict) -> dict:
    """
    Lift over ResultData from 1.1.3 to 1.2.0
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

            # the list of categories which are being treated as support for this run
            variant['var_data']['support_categories'] = config_retrieve(['ValidateMOI', 'support_categories'], [])
            # the list of categories being ignored for this run
            variant['var_data']['ignored_categories'] = config_retrieve(['ValidateMOI', 'ignore_categories'], [])

    data_dict['version'] = '1.2.0'
    return data_dict
