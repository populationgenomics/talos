#!/usr/bin/env python3

"""
A standalone script to read some talos report files, and summarise the number of affected samples involved
"""

from argparse import ArgumentParser
from collections import defaultdict

from cpg_utils import to_path

from talos.models import ResultData

from talos.utils import (
    GeneDict,
    canonical_contigs_from_vcf,
    filter_results,
    find_comp_hets,
    make_flexible_pedigree,
    read_json_from_path,
)


def main(input_json: str, output: str, prefix: bool):
    """
    read the target report, and summarise the number of affected samples involved

    Args:
        input_json ():
        output ():
        prefix ():

    Returns:

    """


def cli_main():
    """
    CLI entrypoint
    """

    parser = ArgumentParser()
    parser.add_argument('--input', help='Path to the input report file')
    parser.add_argument('--output', help='Where to write the output', required=True)
    parser.add_argument('--prefix', action='store_true', help='Expect a consistent Prefix for families')
    args = parser.parse_args()
    main(
        input_json=args.input,
        output=args.output,
        prefix=args.prefix,
    )


if __name__ == '__main__':
    cli_main()
