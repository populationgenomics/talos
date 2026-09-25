"""
code for lifting over models from 2.4.0 to 2.5.0
"""


def resultdata(data_dict: dict) -> dict:
    for res in data_dict['results'].values():
        for var in res['variants']:
            # placeholder to represent that the last run did not record this stat
            var['max_confidence'] = -1
            var['confidence_increase'] = False
    data_dict['version'] = '2.5.0'
    return data_dict


def panelapp(data_dict: dict) -> dict:
    _str_genes = data_dict.pop('str_genes')
    _str_symbols = data_dict.pop('str_symbols')
    data_dict['strs'] = {}
    data_dict['version'] = '2.5.0'
    return data_dict
