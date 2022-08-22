"""
Processes the ACMG additional findings content into a useable config
Should be interchangeable with the panelapp content for downstream

n.b.
the original form of this spreadsheet is horrible to programmatically parse
trim off the top and bottom rows of pure text before parsing
 - i.e. top line should be column headings, ending at last line of gene data
from the original, correct "Disease/Phentyope" to "Disease/Phenotype"
 - corrected form used here, assuming ACMG will correct in future

This format is just not super friendly for parsing, especially with a couple of
per-gene exceptions to the interesting genes. This will be used for MOI tests,
with the Hail layer requiring a messy bit of hard coding of conditionals for now
"""


import json
from csv import DictReader

from argparse import ArgumentParser
from cloudpathlib import AnyPath

KEYS = [
    'Gene',
    'Gene MIM',
    'Disease/Phentyope',
    'Disorder MIM',
    'Phenotype Category',
    'Inheritance',
    'SF List Version',
    'Variants to report',
]
USEFUL_KEYS = {'symbol': 'Gene', 'moi': 'Inheritance'}
FLAG_KEYS = ['Disease/Phentyope', 'Phenotype Category']
MOI_TRANSLATION = {
    'AD': 'Monoallelic',
    'AR': 'Biallelic',
    'XL': 'X-LINKED',
}
LENIENT = 'BOTH'
X_LENIENT = 'X-LINKED'


def main(input_file: str, output_file: str):
    """
    process the input file
    """

    parsed_content = {'metadata': {'panel_name': 'ACMG additional findings'}}
    for row_dict in DictReader(AnyPath(input_file).open()):

        data_blob = {
            key: row_dict[value].rstrip() for key, value in USEFUL_KEYS.items()
        }
        data_blob['flags'] = [row_dict[value].rstrip() for value in FLAG_KEYS]
        data_blob['moi'] = MOI_TRANSLATION[data_blob['moi']]

        if 'only' in row_dict.get('Variants to report'):
            print(
                f'Specific Variants: '
                f'{data_blob["symbol"]} - {row_dict["Variants to report"]}'
            )
        # don't duplicate gene entries
        if data_blob['symbol'] in parsed_content:
            parsed_data = parsed_content[data_blob['symbol']]
            if data_blob['moi'] != parsed_data['moi']:
                if 'x' in parsed_data['moi']:
                    parsed_data['moi'] = X_LENIENT
                else:
                    parsed_data['moi'] = LENIENT
            parsed_data['flags'] = list(
                set(data_blob['flags']).union(set(parsed_data['flags']))
            )

        else:
            parsed_content[data_blob['symbol']] = data_blob

    with AnyPath(output_file).open('w') as handle:
        json.dump(parsed_content, handle, indent=True)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-i', help='input config CSV')
    parser.add_argument('-o', help='output json path')
    args = parser.parse_args()
    main(input_file=args.i, output_file=args.o)
