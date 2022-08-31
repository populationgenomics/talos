"""
script testing methods within reanalysis/validate_categories.py
"""


from reanalysis.validate_categories import update_result_meta


def test_update_results_meta(peddy_ped):
    """
    testing the dict update
    """

    results = {}
    config = {'input_file': 'foo', 'latest_run': 'bar'}
    panelapp = {
        'metadata': {
            'panel_name': 'biff',
            'panel_version': 'pow',
            'additional_panels': [
                {'panel_name': 'extra_panel', 'panel_version': 'extra_version'},
            ],
        }
    }

    big_results = update_result_meta(
        results=results, config=config, pedigree=peddy_ped, panelapp=panelapp
    )
    assert big_results == {
        'metadata': {
            'run_datetime': 'bar',
            'input_file': 'foo',
            'family_breakdown': {
                'affected': 2,
                'male': 3,
                'female': 3,
                'trios': 2,
                '3': 2,
            },
            'panels': [
                {'panel_name': 'extra_panel', 'panel_version': 'extra_version'},
                {'panel_name': 'biff', 'panel_version': 'pow'},
            ],
        }
    }
