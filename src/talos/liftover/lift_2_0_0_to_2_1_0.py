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
            variant['reasons'] = variant['reasons'].pop()
    _categories = data_dict['metadata'].pop('categories')
    data_dict['metadata']['variant_breakdown'] = {}
    data_dict['version'] = '2.1.0'
    return data_dict
