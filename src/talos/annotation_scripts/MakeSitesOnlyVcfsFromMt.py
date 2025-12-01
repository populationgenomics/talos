"""
Write sites-only representations of the original dataset, fragmented as one per MT partition.
"""

from argparse import ArgumentParser

import hail as hl

# VCFs can only be exported from the MatrixTable when all FILTERS entries are accounted for in the header
# Some quality filtering parameters are assigned but not described in the header.
# The solutions are either add an entry in this `filter` dict for each possible FILTER, and add metadata=METADATA in the
# `export_vcf()` method call, or hard filter out quality failures with `mt = mt.filter_rows(mt.filters.length() == 0)`
METADATA = {
    'filter': {
        'LowDepth': {
            'Description': 'LowDepth',
        },
    },
}


def main(input_mt: str, output_dir: str) -> None:
    hl.context.init_spark(master='local[*]', default_reference='GRCh38')
    mt = hl.read_matrix_table(input_mt)

    # remove any filtered out rows
    # todo: maybe remove this later?
    mt = mt.filter_rows(mt.filters.length() == 0)
    hl.export_vcf(mt.rows(), output_dir, parallel='header_per_shard')


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', help='path to input MT file')
    parser.add_argument('--output', help='directory to write VCF fragments to')
    args = parser.parse_args()
    main(input_mt=args.input, output_dir=args.output)
