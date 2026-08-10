"""
code for lifting over models from 2.3.0 to 2.4.0
"""


def panelapp(data_dict: dict) -> dict:
    for gene_data in data_dict['genes'].values():
        gene_data['location'] = ''
    data_dict['version'] = '2.4.0'
    return data_dict


def dl_panelapp(data_dict: dict) -> dict:
    for gene_data in data_dict['genes'].values():
        _mane_symbol = gene_data.pop('mane_symbol', None)
        # cached data only kept the contig, so the full location can't be recovered - re-download for coordinates
        gene_data['location'] = ''
    data_dict['version'] = '2.4.0'
    return data_dict
