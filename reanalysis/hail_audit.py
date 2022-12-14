"""
pulled out the methods relating to the 'hail audit'
mostly just to reduce the line count
partially to make this available in other places
"""


import logging

import hail as hl


def fields_audit(
    mt: hl.MatrixTable, base_fields: list[tuple], nested_fields: dict
) -> bool:
    """
    checks that the required fields are all present before continuing
    """
    problems = []
    # iterate over top-level attributes
    for annotation, datatype in base_fields:
        if annotation in mt.row_value or annotation in mt.row_key:
            if not isinstance(mt[annotation], datatype):
                problems.append(f'{annotation}: {datatype}/{type(mt[annotation])}')
        else:
            problems.append(f'{annotation}:missing')

    for field_group, group_types in nested_fields.items():
        if field_group not in mt.row_value:
            problems.append(f'{field_group}:missing')
        else:
            for annotation, datatype in group_types:
                if annotation in mt[field_group]:
                    if not isinstance(mt[field_group][annotation], datatype):
                        problems.append(
                            f'{annotation}:'
                            f'{datatype}/'
                            f'{type(mt[field_group][annotation])}'
                        )
                else:
                    problems.append(f'{annotation}:missing')
    for problem in problems:
        logging.error(f'MT field: \t{problem}')
    return len(problems) == 0


def vep_audit(mt: hl.MatrixTable, expected_fields: list[tuple]) -> bool:
    """
    check that the required VEP annotations are present
    True if the 'audit' passes (all required fields present)
    """

    problems = []
    # now the content of the transcript_consequences
    if 'vep' not in mt.row_value:
        problems.append('VEP:missing')
    elif 'transcript_consequences' not in mt.vep:
        problems.append('transcript_consequences:missing')
    else:
        fields_and_types = dict(mt.vep.transcript_consequences[0].items())
        for field, field_type in expected_fields:
            if field in fields_and_types:
                if not isinstance(fields_and_types[field], field_type):
                    problems.append(
                        f'{field}:{field_type}/{type(fields_and_types[field])}'
                    )
            else:
                problems.append(f'{field}:missing')

    for problem in problems:
        logging.error(f'VEP field: {problem}')

    return len(problems) == 0
