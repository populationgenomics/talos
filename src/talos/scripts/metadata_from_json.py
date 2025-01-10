#!/usr/bin/env python3

"""
A standalone script to read some talos report files, and summarise the number of affected samples involved

very very simple

Allows for the printing of affected participants grouped by a portion of the family ID
Examples are when the family  ID is prefixed with the year, e.g. 19DNAXXXX for 2019
The simple count function doesn't have this level of granularity

If there is no common prefix, we don't attempt this grouping
"""

import json

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


def main(input_path: str, output_path: str, prefix: int | None = None):
    """
    read the target report, and summarise the number of affected samples involved

    Args:
        input_path (str): where to read the report from
        output_path ():
        prefix (int): optional, if we want to break down by family prefix
    """

    # read the report file, local or cloud
    report = read_json_from_path(input_path, return_model=ResultData)

    # this is a simple overview
    family_breakdown = report.metadata.family_breakdown

    if prefix:
        # setup a section in the dictionary for this
        family_breakdown['grouped_families'] = defaultdict(int)
        for proband in report.results.values():
            family_breakdown[proband.metadata.family_id[:prefix]] += 1

    print(json.dumps(family_breakdown, indent=4))

    # write the output to file
    with to_path(output_path).open('w') as handle:
        json.dump(family_breakdown, handle, indent=4)


def cli_main():
    """
    CLI entrypoint
    """

    parser = ArgumentParser()
    parser.add_argument('--input', help='Path to the input report file')
    parser.add_argument('--output', help='Where to write the output', required=True)
    parser.add_argument('--prefix', type=int, help='breakdown by family prefix', default=None)
    args = parser.parse_args()
    main(
        input_path=args.input,
        output_path=args.output,
        prefix=args.prefix,
    )


if __name__ == '__main__':
    cli_main()
