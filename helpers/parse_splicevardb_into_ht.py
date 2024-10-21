"""
Script for pre-processing the SpliceVarDB tsv file into a Hail table

SpliceVarDB: A comprehensive database of experimentally validated human splicing variants
Sullivan, Patricia J. et al.
The American Journal of Human Genetics, Volume 111, Issue 10, 2164 - 2175

This input file is obtained from the "Download all variants" button on https://compbio.ccia.org.au/splicevardb/

The intention is to use this pre-validated data to accompany and improve upon SpliceAi within Talos.
"""

import gzip
from argparse import ArgumentParser

from talos.utils import get_random_string

from parse_amissense_into_ht import hail_table_from_tsv


"""
Expecting the header fields:
variant_id	hg19	hg38	gene	hgvs	method	classification	location	doi
"""


# important fields
POSITION: str = 'hg38'
GENE: str = 'gene'
CLASSIFICATION: str = 'classification'
METHOD: str = 'method'
VAR_ID: str = 'variant_id'
LOCATION: str = 'location'


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--input', help='input SpliceVarDB tsv file')
    parser.add_argument('--output', help='output SpliceVarDB HT')
    args = parser.parse_args()
    main(args.input, args.output)


def main(input_tsv, output_ht) -> None:
    """
    Main function for SpliceVarDB to HT conversion
    Args:
        input_tsv ():
        output_ht ():
    """

    new_tsv = parse_splicevardb_into_tsv(input_tsv)
    hail_table_from_tsv(new_tsv, output_ht)


def parse_splicevardb_into_tsv(input_tsv: str) -> str:
    """
    Read the SpliceVarDB tsv file and clean it

    Args:
        input_tsv (str): path to the SpliceVarDB tsv file

    Returns:
        file path to the new tsv.gz file
    """

    # holder for new output lines
    new_lines: list[str] = ['chrom\tpos\tref\talt\tgene\tvar_id\tmethod\tclassification\tlocation']
    with open(input_tsv, encoding='utf-8') as handle:
        # file has a header, with no quotes
        header = handle.readline().strip().split('\t')
        for line in handle:
            # there's loads of quotes, get rid of 'em
            stripped_line = line.strip().replace('"', '')
            stripped_list = stripped_line.split('\t')

            # get the position
            position = stripped_list[header.index(POSITION)]
            chrom, pos, ref, alt = position.split('-')

            # update for GRCh38
            if 'chr' not in chrom:
                chrom = f'chr{chrom}'

            # get the gene
            gene = stripped_list[header.index(GENE)]

            # get the variant ID
            variant_id = stripped_list[header.index(VAR_ID)]

            # get the method
            method = stripped_list[header.index(METHOD)]

            # get the location
            location = stripped_list[header.index(LOCATION)]

            # get the location
            classification = stripped_list[header.index(CLASSIFICATION)]
            new_line = f'{chrom}\t{pos}\t{ref}\t{alt}\t{gene}\t{variant_id}\t{method}\t{classification}\t{location}'
            new_lines.append(new_line)

    random_intermediate_file: str = f'{get_random_string()}.tsv.gz'
    with gzip.open(random_intermediate_file, 'wt') as write_handle:
        write_handle.write('\n'.join(new_lines))

    return random_intermediate_file


if __name__ == '__main__':
    cli_main()
