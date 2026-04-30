"""
takes the MitoTip tsv as input, outputs a new VCF file.

input columns:
Position        rCRS    Alt     MitoTIP_Score   Quartile        Count   Percentage      Mitomap_Status
"""

import gzip
from argparse import ArgumentParser
from csv import DictReader
from importlib import resources


def main(input_mitotip: str, output: str):
    with (
        open(input_mitotip, 'r') as handle,
        gzip.open(output, 'wt') as out,
        (resources.files('talos') / 'vcf_headers' / 'mitotip_header.txt').open() as head_in,
    ):
        for line in head_in:
            out.write(line)

        dict_reader = DictReader(handle, delimiter='\t')

        for tsv_line in dict_reader:
            alt = tsv_line['Alt']

            # not sure how to handle this ALT allele, skip these lines
            if alt == ':':
                continue

            # pull out the relevant fields to bake into a VCF
            position = tsv_line['Position']
            ref = tsv_line['rCRS']
            score = tsv_line['MitoTIP_Score']
            mm_status = tsv_line['Mitomap_Status']

            # write a VCF-format row
            out.write(f'chrM\t{position}\t{ref}\t.\t{alt}\t60\tPASS\tmitotip={score};mitomap_status={mm_status}\n')


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', help='Input mitotip text file')
    parser.add_argument('--output', help='output VCF')
    args = parser.parse_args()
    main(input_mitotip=args.input, output=args.output)
