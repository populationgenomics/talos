#!/usr/bin/env python3

"""
Reformat the MANE summary file into a JSON lookup keyed on ENST
https://ftp.ncbi.nlm.nih.gov/refseq/MANE/README.txt
"""

import gzip
import json
from argparse import ArgumentParser
from collections import defaultdict
from csv import DictReader


def mane_to_json(input_path: str, output_path: str):
    """
    Reformat the MANE summary file into a JSON file
    I can't quite grok the annotation from table, so I'm doing the dumb version
    """
    transcript_dict: dict[str, dict[str, str]] = defaultdict(dict)

    with gzip.open(input_path, 'rt') as handle:
        reader = DictReader(handle, delimiter='\t')
        for line in reader:
            enst = line['Ensembl_nuc'].split('.')[0]  # ENST
            ensg = line['Ensembl_Gene'].split('.')[0]  # ENSG
            ensp = line['Ensembl_prot'].split('.')[0]  # ENSP
            symbol = line['symbol'].split('.')[0]  # symbol
            mane_id = line['RefSeq_nuc'].split('.')[0]  # NMID
            mane_status = line['MANE_status']  # "MANE Select" or "MANE Plus Clinical"

            transcript_dict[enst] = {
                'mane_status': mane_status,
                'ensp': ensp,
                'mane_id': mane_id,
                'ensg': ensg,
                'symbol': symbol,
            }

    with open(output_path, 'w') as handle:
        json.dump(transcript_dict, handle, indent=2)


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--input', help='Path to the MANE summary file')
    parser.add_argument('--output', help='Path to write the resulting JSON')
    args = parser.parse_args()
    mane_to_json(input_path=args.input, output_path=args.output)


if __name__ == '__main__':
    cli_main()
