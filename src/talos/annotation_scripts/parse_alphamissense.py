"""
takes an alphamissense tsv as input, filters regions, and outputs a new... VCF file.

input columns:
#CHROM POS REF ALT genome uniprot_id transcript_id protein_variant am_pathogenicity am_class

Using https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz?download=1
- 613MB, containing all the pre-computed data we're interested in
"""

import gzip
from argparse import ArgumentParser
from importlib import resources

resources.files()

def main(input_am: str, output: str):
    with (
        gzip.open(input_am, 'rt') as handle,
        gzip.open(output, 'wt') as out,
        (resources.files('talos') / 'am_header.txt').open() as head_in,
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
    parser.add_argument('--input', help='input gzipped alphamissense tsv')
    parser.add_argument('--output', help='output VCF')
    args = parser.parse_args()
    main(input_am=args.input, output=args.output)
