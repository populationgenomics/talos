"""
code for lifting over models from 2.3.0 to 2.4.0
"""


def resultdata(data_dict: dict) -> dict:
    """
    Lift over ResultData from 2.3.0 to 2.4.0

    The SmallVariant/StructuralVariant union is now discriminated on an explicit `kind` tag
    rather than on the presence of `transcript_consequences`. Historical documents predate
    the tag, so inject it using the same signal the old implicit discrimination relied on:
    a variant carrying `transcript_consequences` is a small variant, otherwise structural.
    """
    for result in data_dict['results'].values():
        for variant in result['variants']:
            var_data = variant['var_data']
            var_data['kind'] = 'small' if 'transcript_consequences' in var_data else 'sv'

    data_dict['version'] = '2.4.0'
    return data_dict
