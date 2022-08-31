"""
script testing methods within reanalysis/validate_categories.py
"""


from reanalysis.validate_categories import update_result_meta


def test_update_results_meta(peddy_ped):
    """
    testing the dict update
    """

    results = {}
    config = {'input_file': 'foo', 'latest_run': 'bar', 'cohort': 'cohort'}
    panelapp = {
        'metadata': {
            'panel_name': 'biff',
            'panel_version': 'pow',
            'panel_id': 'wallop',
            'additional_panels': [
                {
                    'panel_name': 'extra_panel',
                    'panel_version': 'extra_version',
                    'panel_id': 2,
                },
            ],
        }
    }
    ped_samples = ['male', 'female', 'mother_1', 'father_1', 'mother_2', 'father_2']

    big_results = update_result_meta(
        results=results,
        config=config,
        pedigree=peddy_ped,
        panelapp=panelapp,
        samples=ped_samples,
    )
    assert big_results == {
        'metadata': {
            'run_datetime': 'bar',
            'input_file': 'foo',
            'cohort': 'cohort',
            'family_breakdown': {
                'affected': 2,
                'male': 3,
                'female': 3,
                'trios': 2,
                '3': 2,
            },
            'panels': [
                {
                    'panel_name': 'extra_panel',
                    'panel_version': 'extra_version',
                    'panel_id': 2,
                },
                {'panel_name': 'biff', 'panel_version': 'pow', 'panel_id': 'wallop'},
            ],
        }
    }
