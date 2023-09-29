"""
Takes an input MT, and extracts a VCF-format representation.

This is currently required as the end-to-end CPG pipeline doesn't currently
store intermediate files. To simulate workflows running on VCF files, we
have to regenerate a VCF representation from a MT.

Hard coded additional header file for VQSR content
When Hail extracts a VCF from a MT, it doesn't contain any custom field
definitions, e.g. 'VQSR' as a Filter field. This argument allows us to
specify additional lines which are required to make the final output valid
within the VCF specification

If the file was not processed with VQSR, there are no negatives to including
this additional header line
"""


from argparse import ArgumentParser

import hail as hl

from cpg_utils import to_path
from cpg_utils.hail_batch import init_batch, output_path


def main(mt_path: str, write_path: str):
    """
    takes an input MT, and reads it out as a VCF
    inserted new conditions to minimise the data produced

    Args:
        mt_path ():
        write_path ():
    """
    init_batch()

    mt = hl.read_matrix_table(mt_path)

    # this temp file needs to be in GCP, not local
    # otherwise the batch that generates the file won't be able to read
    additional_cloud_path = output_path('additional_header.txt', 'tmp')

    with to_path(additional_cloud_path).open('w') as handle:
        handle.write('##FILTER=<ID=VQSR,Description="VQSR triggered">')

    # remove potentially problematic field from gVCF
    if 'gvcf_info' in mt.row_value:
        mt = mt.drop('gvcf_info')

    # filter out filter failures and non-variant rows (prior to VEP)
    mt = mt.filter_rows(mt.filters.length() == 0)
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.n_non_ref > 0)

    # required to repeat the MT -> VCF -> MT cycle
    if 'info' in mt.row_value and 'AC' in mt.info:
        if isinstance(mt['info']['AC'], hl.Int32Expression):
            mt = mt.annotate_rows(info=mt.info.annotate(AC=[mt.info.AC]))
        else:
            mt = mt.annotate_rows(info=mt.info.annotate(AC=[1]))

    hl.export_vcf(
        mt,
        write_path,
        append_to_header=additional_cloud_path,
        tabix=True,
    )


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', type=str, help='input MatrixTable path')
    parser.add_argument('--output', type=str, help='path to write VCF out to')
    args = parser.parse_args()
    main(mt_path=args.input, write_path=args.output)
