"""
code for lifting over models from 2.1.0 to 2.2.0
"""


def panelapp(data_dict: dict) -> dict:
    data_dict |= {'str_genes': set(), 'str_symbols': set()}
    data_dict['version'] = '2.3.0'
    return data_dict


def resultdata(data_dict: dict) -> dict:
    for result in data_dict['results'].values():
        for variant in result['variants']:
            variant['categories'] = dict.fromkeys(variant['categories'], variant['evidence_last_updated'])
            _indi = variant.pop('independent')
    data_dict['version'] = '2.2.0'
    return data_dict
