"""
code for lifting over models from None/unspecified to 1.0.0
"""


def lift_phenotype_matched_panels(data_dict: dict) -> dict:
    """
    Lift over PhenotypeMatchedPanels from None to 1.0.0
    requires the migration of HPO terms from strings to dicts
    """
    return data_dict


def lift_panel_app(data_dict: dict) -> dict:
    """
    Lift over PanelApp from None to 1.0.0
    """
    return data_dict


def lift_result_data(data_dict: dict) -> dict:
    """
    Lift over PanelApp from None to 1.0.0
    """
    return data_dict


def lift_historic_panels(data_dict: dict) -> dict:
    """
    Lift over HistoricPanels from None to 1.0.0
    """
    return data_dict


def lift_historic_variants(data_dict: dict) -> dict:
    """
    Lift over HistoricVariants from None to 1.0.0
    """
    return data_dict
