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

vars = {
    'chr1:21706892',
    'chr2:135912503',
    'chr2:135920591',
    'chr6:26090951',
    'chr6:39887558',
    'chr6:52043102',
    'chr6:52043699',
    'chr11:32392032',
    'chr11:62841775',
    'chr11:123057736',
    'chr12:89470359',
    'chr12:120291834',
    'chr16:89279566',
    'chrX:71109321',
}


def main(input_am: str, output: str):
    with (
        gzip.open(input_am, 'rt') as handle,
        gzip.open(output, 'wt') as out,
        (resources.files('talos') / 'vcf_headers' / 'am_header.txt').open() as head_in,
    ):
        for line in head_in:
            out.write(line)

        for line in handle:
            if line.startswith('#'):
                continue

            llist = line.rstrip().split()

            if f'{llist[0]}:{llist[1]}' not in vars:
                continue

            out.write(
                f'{llist[0]}\t{llist[1]}\t.\t{llist[2]}\t{llist[3]}\t60\tPASS\tam_class={llist[9]};am_score={llist[8]};am_transcript={llist[6].split(".")[0]}\n',
            )


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', help='input gzipped alphamissense tsv')
    parser.add_argument('--output', help='output VCF')
    args = parser.parse_args()
    main(input_am=args.input, output=args.output)
