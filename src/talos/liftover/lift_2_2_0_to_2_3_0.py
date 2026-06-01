"""
code for lifting over models from 2.2.0 to 2.3.0
"""


def dl_panelapp(data_dict: dict) -> dict:
    for _, gene_data in data_dict['genes'].items():
        gene_data['panel_confidences'] = {}
        for _, panel_data in gene_data['panels'].items():
            panel_data['confidence'] = 3
    data_dict['version'] = '2.3.0'
    return data_dict
