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
from collections import defaultdict

import toml


def coord_to_string(coord: dict) -> str:
    """
    converts a coordinate dict to a string

    Args:
        coord (dict): a coordinate dict

    Returns:
        str: a string representation of the coordinate dict
    """
    return f"{coord['chrom']}-{coord['pos']}-{coord['ref']}-{coord['alt']}"


def main(input_file: str, output: str, config: str):
    """
    reads in the input file, shrinks it, and writes the output file

    Args:
        input_file (str):
        output (str):
        config (str):

    Returns:

    """
    with open(config, encoding='utf-8') as f:
        config = toml.load(f)

    with open(input_file, encoding='utf-8') as f:
        data = json.load(f)

    lil_data = {
        'metadata': {'categories': config['categories']},
        'results': defaultdict(dict),
    }

    for individual, details in data['results'].items():
        for variant in details['variants']:
            var_data = variant['var_data']
            lil_data['results'][individual][coord_to_string(var_data['coords'])] = {
                'categories': var_data['categories'],
                # 'labels': variant['labels'],
                'support_vars': variant['support_vars'],
                # 'independent': variant['independent'],
            }

    with open(output, 'w', encoding='utf-8') as f:
        json.dump(lil_data, f, indent=4)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('input_file', help='the input file to process')
    parser.add_argument('output_file', help='the output file to write to')
    parser.add_argument('config_file', help='the config file to use')
    args = parser.parse_args()

    main(input_file=args.input_file, output=args.output_file, config=args.config_file)
