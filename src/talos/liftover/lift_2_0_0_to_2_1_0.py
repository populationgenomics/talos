"""
code for lifting over models from 1.2.0 to 2.0.0
"""


def panelapp(data_dict: dict) -> dict:
    for participant_data in data_dict['participants'].values():
        _ext = participant_data.pop('ext_id')

    data_dict['version'] = '2.1.0'
    return data_dict


def resultdata(data_dict: dict) -> dict:
    print(data_dict)
    for result in data_dict['results'].values():
        for variant in result['variants']:
            variant['reasons'] = variant['reasons'].pop()
    _categories = data_dict['metadata'].pop('categories')
    data_dict['metadata']['variant_breakdown'] = {}
    data_dict['version'] = '2.1.0'
    return data_dict


def historicvariants(data_dict: dict) -> dict:
    """
    Lift over HistoricVariants from 2.0.0 to 2.1.0
    requires the translation of the category short IDs (e.g. '1', '3', '4') to the full names
    """

    # this was the exact lookup present as-of 2.1.0
    category_translator: dict[str, str] = {
        '1': 'ClinVarP/LP',
        '3': 'HighImpact',
        '4': 'DeNovo',
        '5': 'SpliceAI',
        '6': 'AlphaMissense',
        'pm5': 'PM5',
        'sv1': 'LofSv',
        'svdb': 'SpliceVarDB',
        'exomiser': 'Exomiser',
    }

    for results in data_dict['results'].values():
        for details in results.values():
            details['categories'] = {
                category_translator.get(cat_id, cat_id): date for cat_id, date in details['categories'].items()
            }

    _cat_dict = data_dict.pop('metadata')

    data_dict['version'] = '2.1.0'
    return data_dict
