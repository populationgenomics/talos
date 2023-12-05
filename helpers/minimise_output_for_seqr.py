"""
a script to post-process the AIP output data

produces a minimised representation of the output data,
containing only the data required for the SEQR app.

- Metadata, containing the description for each flag
- Results, containing per-individual results
    - Individual ID
        - Variant ID
            - Categories (list)
            - Labels (list)
            - Support Variants (list)
            - Independent (bool)
"""

import json
from argparse import ArgumentParser

from reanalysis.models import MiniForSeqr, MiniVariant, ResultData


def coord_to_string(coord: dict) -> str:
    """
    converts a coordinate dict to a string

    Args:
        coord (dict): a coordinate dict

    Returns:
        str: a string representation of the coordinate dict
    """
    return f"{coord['chrom']}-{coord['pos']}-{coord['ref']}-{coord['alt']}"


def main(input_file: str, output: str):
    """
    reads in the input file, shrinks it, and writes the output file

    Args:
        input_file (str):
        output (str):
    """

    with open(input_file, encoding='utf-8') as f:
        data = ResultData.model_validate(json.load(f))

    lil_data = MiniForSeqr(
        **{
            'metadata': {'categories': data.metadata.categories},
        }
    )

    for individual, details in data.results.items():
        lil_data.results[individual] = {}
        for variant in details.variants:
            var_data = variant.var_data
            lil_data.results[individual][
                var_data.coordinates.string_format
            ] = MiniVariant(
                **{
                    'categories': variant.categories,
                    'support_vars': variant.support_vars,
                    'independent': variant.independent,
                }
            )

    with open(output, 'w', encoding='utf-8') as f:
        f.write(MiniForSeqr.model_validate(lil_data).model_dump_json(indent=4))


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('input_file', help='the input file to process')
    parser.add_argument('output_file', help='the output file to write to')
    args = parser.parse_args()

    main(input_file=args.input_file, output=args.output_file)
