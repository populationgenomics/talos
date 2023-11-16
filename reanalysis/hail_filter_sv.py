"""
A Hail filtering process for labelling analysis-relevant SVs
Initially this will only contain a single category

CategoryBooleanSV1:
- rare
- deletion in a listed gene
- ???
"""


from argparse import ArgumentParser

import hail as hl

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import init_batch, genome_build

from reanalysis.hail_filter_and_label import (
    green_and_new_from_panelapp,
    subselect_mt_to_pedigree,
)
from reanalysis.utils import get_logger, read_json_from_path


def main(
    mt_path: str,
    panelapp_path: str,
    pedigree: str,
    vcf_out: str,
    dataset: str | None = None,
):
    """
    Read MT, filter, and apply category annotation
    Export as a VCF

    Args:
        mt_path (str): where to find vcf output
        panelapp_path ():
        pedigree ():
        vcf_out (str): where to write VCF out
        dataset (str): optional dataset name to write output for
    """
    dataset = dataset or get_config(True)['workflow']['dataset']

    # this workflow should be too simple (at this time) for checkpoints

    # read the parsed panelapp data
    get_logger().info(f'Reading PanelApp data from {panelapp_path!r}')
    panelapp = read_json_from_path(panelapp_path)['genes']  # type: ignore

    # pull green and new genes from the panelapp data
    green_expression, new_expression = green_and_new_from_panelapp(panelapp)

    # initiate Hail with defined driver spec.
    get_logger().info(f'Starting Hail with reference genome {genome_build()}')
    init_batch(driver_cores=8, driver_memory='highmem')

    # if we already generated the annotated output, load instead
    if not to_path(mt_path.rstrip('/') + '/').exists():
        raise FileExistsError(f'Input MatrixTable doesn\'t exist: {mt_path}')

    mt = hl.read_matrix_table(mt_path)

    # subset to currently considered samples
    mt = subselect_mt_to_pedigree(mt, pedigree=pedigree)


if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('--mt', required=True, help='path to input MT')
    parser.add_argument('--panelapp', type=str, required=True, help='panelapp JSON')
    parser.add_argument('--pedigree', type=str, required=True, help='Cohort Pedigree')
    parser.add_argument('--vcf_out', help='Where to write the VCF', required=True)
    parser.add_argument('--dataset', help='Dataset to write output for')
    args = parser.parse_args()

    logger = get_logger(__name__)
    logger.info(f'Running Hail filtering process for {args.dataset} SVs')

    main(
        mt_path=args.mt,
        panelapp_path=args.panelapp,
        pedigree=args.pedigree,
        vcf_out=args.vcf_out,
        dataset=args.dataset,
    )
