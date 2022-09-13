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


from csv import DictReader
import json
import logging
import os
import re
import sys

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


def get_data_from_row(row_data: dict) -> dict:
    """

    Parameters
    ----------
    row_data : the DictReader entry for the row

    Returns
    -------
    a dictionary of that content, parsed
    """

    # parse all top level data
    data_blob = {key: row_data[value].rstrip() for key, value in USEFUL_KEYS.items()}

    # assign flags to the row as appropriate
    data_blob['flags'] = [row_data[value].rstrip() for value in FLAG_KEYS]

    # replace the default MOI with a 'simple' MOI
    data_blob['moi'] = MOI_TRANSLATION[data_blob['moi']]

    # find any specific targets
    data_blob['specific_type'] = re.findall(
        VARIANT_TYPE, row_data.get('Variants to report', '')
    )
    data_blob['specific_variant'] = re.findall(
        SPECIFIC_CHANGE, row_data.get('Variants to report', '')
    )

    if 'only' in row_data.get('Variants to report', ''):
        print(
            f'Specific Variants: '
            f'{data_blob["symbol"]} - {row_data["Variants to report"]}'
        )

    return data_blob


def main(input_file: str, output_file: str):
    """
    process the input file
    """

    parsed_content: dict[str, list[dict] | dict] = {
        'metadata': [{'name': 'ACMG additional findings', 'version': 1.0}]
    }
    gene_id_map = read_json_from_path(ADDITIONAL_MAP)
    for row_dict in DictReader(to_path(input_file).open()):

        data_blob = get_data_from_row(row_dict)

        # don't duplicate gene entries - use ESG as key
        ensg_id = gene_id_map[data_blob['symbol']]

        # if we already found an entry for this gene, squash them
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
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )
    parser = ArgumentParser()
    parser.add_argument('-i', help='input config CSV', required=True)
    parser.add_argument('-o', help='output json path', required=True)
    args = parser.parse_args()
    main(input_file=args.i, output_file=args.o)
