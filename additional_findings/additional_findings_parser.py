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

The idea here will be to have a user step in where specific variants or types
are identified. The two types will be
"""


from csv import DictReader
import json
import os
import re

from argparse import ArgumentParser

from cpg_utils import to_path
from reanalysis.utils import read_json_from_path


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
VARIANT_TYPE = re.compile(r'\((\w+) variants only\)')
SPECIFIC_CHANGE = re.compile(r'(p.[A-Z][0-9]+[A-Z]) (?:\w+ )?only')

# placeholder, find a better way to get in future
ADDITIONAL_MAP = os.path.join(os.path.dirname(__file__), 'additional_map.json')


def main(input_file: str, output_file: str):
    """
    process the input file
    """

    parsed_content: dict[str, list[dict] | dict] = {
        'metadata': [{'name': 'ACMG additional findings', 'version': 1.0}]
    }
    gene_id_map = read_json_from_path(ADDITIONAL_MAP)
    for row_dict in DictReader(to_path(input_file).open()):

        data_blob = {
            key: row_dict[value].rstrip() for key, value in USEFUL_KEYS.items()
        }
        data_blob['flags'] = [row_dict[value].rstrip() for value in FLAG_KEYS]
        data_blob['moi'] = MOI_TRANSLATION[data_blob['moi']]
        data_blob['specific_type'] = re.findall(
            VARIANT_TYPE, row_dict.get('Variants to report', '')
        )
        data_blob['specific_variant'] = re.findall(
            SPECIFIC_CHANGE, row_dict.get('Variants to report', '')
        )

        # don't duplicate gene entries - use ESG as key
        ensg_id = gene_id_map[data_blob['symbol']]

        if 'only' in row_dict.get('Variants to report'):
            print(
                f'Specific Variants: '
                f'{data_blob["symbol"]} ({ensg_id }) - {row_dict["Variants to report"]}'
            )

        if data_blob['symbol'] in parsed_content:
            parsed_data = parsed_content[ensg_id]

            # e.g. if for some reason there are two entries, one AR one AD
            # change the discovered MOI to BOTH
            if data_blob['moi'] != parsed_data['moi']:
                if 'x' in parsed_data['moi']:
                    parsed_data['moi'] = X_LENIENT
                else:
                    parsed_data['moi'] = LENIENT
            parsed_data['flags'] = list(
                set(data_blob['flags']).union(set(parsed_data['flags']))
            )
            parsed_data['specific_variant'].extend(data_blob['specific_variant'])
            parsed_data['specific_type'].extend(data_blob['specific_type'])

        else:
            parsed_content[ensg_id] = data_blob

    with to_path(output_file).open('w') as handle:
        json.dump(parsed_content, handle, indent=True)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-i', help='input config CSV', required=True)
    parser.add_argument('-o', help='output json path', required=True)
    args = parser.parse_args()
    main(input_file=args.i, output_file=args.o)
