"""
takes an input file, works out the format, and
"""


from argparse import ArgumentParser
import hail as hl
from cpg_utils.hail_batch import init_batch


def main(input_mt: str, output_path: str):
    """
    takes an input MT, and reads it out as a VCF
    :param input_mt:
    :param output_path:
    :return:
    """
    init_batch()

    matrix = hl.read_matrix_table(input_mt)

    hl.export_vcf(
        matrix,
        output_path,
        tabix=True,
    )


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument(
        '--input',
        type=str,
        help='input MatrixTable path',
    )
    parser.add_argument('--output', type=str, help='path to write VCF out to')
    args = parser.parse_args()
    main(input_mt=args.input, output_path=args.output)
