"""
Parses a GFF3 file, and generates a BED file of gene regions, plus padding
"""

import gzip
import re
from argparse import ArgumentParser


CHROM_INDEX = 0
RESOURCE_INDEX = 1
TYPE_INDEX = 2
START_INDEX = 3
END_INDEX = 4
DETAILS_INDEX = 8

# +/- this is added to each gene region
FLANKING_REGION = 2000
GENE_ID_RE = re.compile(r'ID=gene:(ENSG\d+);')
GENE_NAME_RE = re.compile(r'Name=([\w-]+);')


def main(gff3_file: str, output: str, flanking: int = 2000):
    """
    Read the GFF3 file, and generate a BED file of gene regions, plus padding
    Args:
        gff3_file (str): path to the GFF3 file
        output (str): path to the intended BED output file
        flanking (int): number of bases to add before and after each gene
    """

    # open and iterate over the GFF3 file
    with gzip.open(gff3_file, 'rt') as handle, open(output, 'w') as write_handle:
        for line in handle:
            # skip over headers and dividing lines
            if line.startswith('#'):
                continue

            line_as_list = line.rstrip().split('\t')

            # skip over non-genes (e.g. pseudogenes, ncRNA)
            # only focus on Ensembl genes/transcripts
            if line_as_list[TYPE_INDEX] != 'gene' or 'ensembl' not in line_as_list[RESOURCE_INDEX]:
                continue

            # extract the gene name from the details field
            # allowing for some situations that don't work,
            # e.g. ENSG00000225931 (novel transcript, to be experimentally confirmed)
            # search for ID and transcript separately, ordering not guaranteed
            try:
                gene_name = GENE_NAME_RE.search(line_as_list[DETAILS_INDEX]).group(1)
                gene_id = GENE_ID_RE.search(line_as_list[DETAILS_INDEX]).group(1)
            except AttributeError:
                print(f'Failed to extract gene name from {line_as_list[DETAILS_INDEX]}')
                continue

            # write the line to the output
            output_list = [
                f'chr{line_as_list[CHROM_INDEX]}',
                str(int(line_as_list[START_INDEX]) - flanking),
                str(int(line_as_list[END_INDEX]) + flanking),
                f'{gene_id};{gene_name}',
            ]
            write_handle.write('\t'.join(output_list) + '\n')


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--gff3', help='Path to the compressed GFF3 file')
    parser.add_argument('--output', help='Path to output file')
    parser.add_argument('--flanking', help='Regions to add to each gene', default=2000)
    args = parser.parse_args()
    main(gff3_file=args.gff3, output=args.output, flanking=args.flanking)


if __name__ == '__main__':
    cli_main()
