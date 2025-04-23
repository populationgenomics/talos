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
