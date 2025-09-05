"""
code for lifting over models from 2.1.0 to 2.2.0
"""


def resultdata(data_dict: dict) -> dict:
    for result in data_dict['results'].values():
        for variant in result['variants']:
            variant['categories'] = {cat: variant['evidence_last_updated'] for cat in variant['categories']}
            _indi = variant.pop('independent')
    data_dict['version'] = '2.2.0'
    return data_dict
