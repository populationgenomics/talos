"""
code for lifting over models from 1.2.0 to 2.0.0
"""


def panelapp(data_dict: dict) -> dict:
    for participant_data in data_dict['participants'].values():
        _ext = participant_data.pop('ext_id')

    data_dict['version'] = '2.1.0'
    return data_dict
