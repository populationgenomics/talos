"""
code for lifting over models from 1.2.0 to 2.0.0
"""


def panelapp(data_dict: dict) -> dict:
    """
    Lift is not valid for PanelApp between 1.2.0 and 2.0.0
    These two objects are fundamentally different between these versions, and should be regenerated from scratch
    """
    _ = data_dict
    raise ValueError('PanelApp data cannot be lifted over between 1.2.0 and 2.0.0 - regenerate the data from scratch')


def resultdata(data_dict: dict) -> dict:
    """
    only substantial change here is the removal of the projects key from the metadata
    this was explicitly tied to seqr, which is no longer a core assumption
    """

    # stop using seqr project(s) as a core value in the data
    _projects = data_dict['metadata'].pop('projects')

    # remove the projects key from the metadata
    data_dict['metadata']['panels'] = {panel['id']: panel for panel in data_dict['metadata']['panels']}

    data_dict['version'] = '2.0.0'
    return data_dict
