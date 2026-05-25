"""
code for lifting over models from 2.2.0 to 2.3.0
"""


def panelapp(data_dict: dict) -> dict:
    data_dict |= {'str_genes': set(), 'str_symbols': set()}
    data_dict['version'] = '2.3.0'
    return data_dict
