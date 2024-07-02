"""
code for lifting over models from 1.0.2 to 1.0.3
"""


def resultdata(data_dict: dict) -> dict:
    """
    Lift over ResultData from 1.0.0 to 1.0.1
    Requires the adjustment of ResultMeta's "container" to "version"
    not using this variable at the moment
    """
    # check we're upgrading the right version
    # this could be from any prior version
    assert data_dict['version'] < '1.0.3'
    metadata = data_dict['metadata']
    assert 'version' not in metadata, metadata
    assert 'container' in metadata, metadata
    metadata['version'] = metadata.pop('container')
    return data_dict


def historicvariants(data_dict: dict) -> dict:
    """
    Lift over HistoricVariants from 1.0.2 to 1.0.3
    Logs the number of ClinVar stars, so that we can identify upgraded stars in future
    without subdividing the Category 1 into multiple categories per star rating
    """
    assert data_dict['version'] < '1.0.3'
    ...
    return data_dict
