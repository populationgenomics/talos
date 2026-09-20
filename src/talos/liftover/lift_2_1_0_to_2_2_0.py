"""
code for lifting over models from 2.1.0 to 2.2.0
"""

from loguru import logger

from talos.static_values import get_granular_date


def dl_panelapp(data_dict: dict) -> dict:
    logger.warning(
        'This PanelApp data does not contain a download date. Talos will warn when PanelApp data is older than 2 '
        'months post-download, but this check will be skipped. We advise you delete this cached data file, and '
        'generate another.',
    )

    data_dict['date'] = get_granular_date()
    data_dict['version'] = '2.2.0'
    return data_dict


def resultdata(data_dict: dict) -> dict:
    for result in data_dict['results'].values():
        for variant in result['variants']:
            vd = variant['var_data']
            categories = variant['categories'] if 'categories' in variant else vd['categories']
            variant['categories'] = dict.fromkeys(categories, variant['evidence_last_updated'])
            if 'independent' in variant:
                _indi = variant.pop('independent')
    data_dict['version'] = '2.2.0'
    return data_dict
