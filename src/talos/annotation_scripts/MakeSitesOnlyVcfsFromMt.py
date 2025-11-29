"""
Write sites-only representations of the original dataset, fragmented as one per MT partition.
"""

from argparse import ArgumentParser

import hail as hl


def main(input_mt: str, output_dir: str) -> None:
    hl.context.init_spark(master='local[*]', default_reference='GRCh38', quiet=True)
    mt = hl.read_matrix_table(input_mt)
    hl.export_vcf(mt.rows(), output_dir, parallel='header_per_shard')


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', help='path to input MT file')
    parser.add_argument('--output', help='directory to write VCF fragments to')
    args = parser.parse_args()
    main(input_mt=args.input, output_dir=args.output)
