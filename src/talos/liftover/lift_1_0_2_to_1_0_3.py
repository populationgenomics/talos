"""
code for lifting over models from 1.0.2 to 1.0.3
"""


def resultdata(data_dict: dict) -> dict:
    """
    Lift over ResultData from 1.0.2 to 1.0.3
    Requires the adjustment of ResultMeta's "container" to "version"
    not using this variable at the moment
    """
    # check we're upgrading the right version
    # this could be from any prior version
    if not data_dict['version'] < '1.0.3':
        raise AssertionError(f'This method cannot upgrade from {data_dict["version"]}')
    metadata = data_dict['metadata']
    print(metadata)

    if ('version' in metadata) or ('container' not in metadata):
        raise AssertionError(f'Metadata should have container, not version: {metadata}')
    metadata['version'] = metadata.pop('container')
    return data_dict


def historicvariants(data_dict: dict) -> dict:
    """
    Lift over HistoricVariants from 1.0.2 to 1.0.3
    Logs the number of ClinVar stars, so that we can identify upgraded stars in future
    without subdividing the Category 1 into multiple categories per star rating
    """
    if not data_dict['version'] < '1.0.3':
        raise AssertionError(f'This method cannot upgrade from {data_dict["version"]}')
    return data_dict
