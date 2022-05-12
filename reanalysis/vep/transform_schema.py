"""
extracting out method for altering schema, for future use
"""


def transform_schema(json_schema: str):
    """
    example input format:
        'Struct{assembly_name:String,allele_string:String...'
    example output format:
        'struct{assembly_name:str,allele_string:str...'
    :param json_schema:
    :return:
    """
    for k, v in [
        ('String', 'str'),
        ('Array', 'array'),
        ('Set', 'set'),
        ('Tuple', 'tuple'),
        ('Struct', 'struct'),
        ('[', '<'),
        (']', '>'),
        ('Int32', 'int32'),
        ('Int64', 'int64'),
        ('Float64', 'float64'),
        ('Float32', 'float32'),
    ]:
        json_schema = json_schema.replace(k, v)
    return json_schema
