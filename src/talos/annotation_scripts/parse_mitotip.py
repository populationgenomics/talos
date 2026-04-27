"""
takes the MitoTip tsv as input, outputs a new VCF file.

input columns:
Position        rCRS    Alt     MitoTIP_Score   Quartile        Count   Percentage      Mitomap_Status
"""

import gzip
from argparse import ArgumentParser
from importlib import resources


def main(input_mitotip: str, output: str):
    with (
        open(input_mitotip, 'r') as handle,
        gzip.open(output, 'wt') as out,
        (resources.files('talos') / 'vcf_headers' / 'mitotip_header.txt').open() as head_in,
    ):
        for line in head_in:
            out.write(line)

        for line in handle:
            if line.startswith('#'):
                continue

            llist = line.rstrip().split()

            out.write(
                f'{llist[0]}\t{llist[1]}\t.\t{llist[2]}\t{llist[3]}\t60\tPASS\tam_class={llist[9]};am_score={llist[8]};am_transcript={llist[6].split(".")[0]}\n',
            )


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', help='Input mitotip text file')
    parser.add_argument('--output', help='output VCF')
    args = parser.parse_args()
    main(input_mitotip=args.input, output=args.output)
