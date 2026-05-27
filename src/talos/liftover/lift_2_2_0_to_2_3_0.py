"""
code for lifting over models from 2.2.0 to 2.3.0
"""

from loguru import logger

from talos.static_values import get_granular_date


def dl_panelapp(data_dict: dict) -> dict:
    for _, gene_data in data_dict['genes'].items():
        for _, panel_data in gene_data['panels'].items():
            panel_data['confidence'] = 3
    data_dict['version'] = '2.3.0'
    return data_dict
