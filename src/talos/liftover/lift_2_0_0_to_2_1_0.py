"""
code for lifting over models from 2.0.0 to 2.1.0
"""


def panelapp(data_dict: dict) -> dict:
    for participant_data in data_dict['participants'].values():
        _ext = participant_data.pop('ext_id')

    data_dict['version'] = '2.1.0'
    return data_dict


def resultdata(data_dict: dict) -> dict:
    for result in data_dict['results'].values():
        for variant in result['variants']:
            # collapse the list of reasons to a single string; an empty list becomes '' (the current default)
            reasons = variant['reasons']
            variant['reasons'] = reasons[-1] if reasons else ''
    _categories = data_dict['metadata'].pop('categories')
    data_dict['metadata']['variant_breakdown'] = {}
    data_dict['version'] = '2.1.0'
    return data_dict
