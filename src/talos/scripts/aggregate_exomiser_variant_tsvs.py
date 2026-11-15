"""
script for taking the per-proband Variant TSVs from Exomiser and aggregating them into a single JSON file.

This was a Talos 1.0 legacy implementation, the script is salvaged here in case non-Hail Exomiser compatibility is
sought in future.

This in no way reimplements or runs Exomiser, and assumes you have already run Exomiser on your data.
The expected data structure and excellent performance has been verified on version 14.0.0

Please see https://github.com/exomiser/Exomiser for more information
"""

import json
from argparse import ArgumentParser
from collections import defaultdict
from csv import DictReader

from cloudpathlib.anypath import to_anypath


ORDERED_ALLELES: list[str] = [f'chr{x}' for x in list(range(1, 23))] + ['chrX', 'chrY', 'chrM']
PROBAND_DICT = dict[str, dict[str, list[dict]]]


def process_tsv(tsv_path: str) -> PROBAND_DICT:
    """
    Read a single TSV file, and return the parsed results as a dict of variant IDs to lists of sample IDs & details

    Args:
        tsv_path (str): path to this TSV file

    Returns:
        a dictionary of variant keys to lists of the 'proband_rank_moi'
    """

    # this is on the assumption that we retain the path/to/file/PROBAND.variant.tsv naming convention
    # if we read into the batch as /tmp/path/PROBAND this logic is still valid
    file_as_path = to_anypath(tsv_path)
    proband = file_as_path.name.split('.')[0]

    proband_dictionary: PROBAND_DICT = {proband: defaultdict(list)}

    with file_as_path.open() as handle:
        reader = DictReader(handle, delimiter='\t')
        for row in reader:
            assert isinstance(row, dict)

            # Exomiser contains "MT" on all genome builds, which Hail does not accept. Overrule this behaviour.
            contig = row['CONTIG']
            if contig == 'MT':
                contig = 'M'

            # this isn't prefixed with 'chr' in the TSV, so we add it
            chrom = f'chr{contig}'
            pos = row['START']
            ref = row['REF']
            alt = row['ALT']

            variant_key = f'{chrom}:{pos}:{ref}:{alt}'

            rank = int(row['#RANK'])
            moi = row['MOI']

            # adds easily parsed details to the family dictionary, easily extensible
            proband_dictionary[proband][variant_key].append(
                {
                    'rank': rank,
                    'moi': moi,
                },
            )

    return proband_dictionary


def main(input_tsvs: list[str], output_path: str):
    """
    Combine the per-proband TSVs into a single JSON file, and write as a Hail Table

    The data will be aggregated per variant, not per proband

    Args:
        input_tsvs (list[str]): list of paths to the per-proband TSVs
        output_path (str): where to write the output
    """

    if not output_path.endswith('.json'):
        raise ValueError('Output path must end with .json')

    all_proband_dictionary: PROBAND_DICT = {}
    for tsv in input_tsvs:
        all_proband_dictionary.update(process_tsv(tsv))

    # write the JSON
    with open(f'{output_path}', 'w') as handle:
        json.dump(all_proband_dictionary, handle, indent=4)


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--input', help='Path to the Variant result TSVs', nargs='+', required=True)
    parser.add_argument('--output', help='Where to write the output json', required=True)
    args = parser.parse_args()

    main(input_tsvs=args.input, output_path=args.output)


if __name__ == '__main__':
    cli_main()
