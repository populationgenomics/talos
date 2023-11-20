"""
A Hail filtering process for labelling analysis-relevant SVs
Initially this will only contain a single category

CategoryBooleanSV1:
- rare
- predicted LoF in a listed gene
"""


from argparse import ArgumentParser

import hail as hl

from cpg_utils import to_path
from cpg_utils.hail_batch import init_batch, genome_build

from reanalysis.hail_filter_and_label import (
    green_and_new_from_panelapp,
    subselect_mt_to_pedigree,
    ONE_INT,
    MISSING_INT,
)
from reanalysis.utils import get_logger, read_json_from_path


def filter_matrix_by_af(
    mt: hl.MatrixTable, af_threshold: float = 0.05
) -> hl.MatrixTable:
    """
    Filter a MatrixTable on AF, allow AF to be missing

    Args:
        mt (hl.MatrixTable): the input MT
        af_threshold (float): filtering threshold in gnomad v2.1

    Returns:
        same MT, with common variants removed
    """

    return mt.filter_rows(
        (hl.or_else(mt.info['gnomad_v2.1_sv_AF'], MISSING_INT) < af_threshold)
    )


def restructure_mt_by_gene(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    split each transcript consequence annotation onto a separate row

    Args:
        mt (hl.MatrixTable): the input MT

    Returns:
        hl.MatrixTable: the filtered MT
    """

    # split out consequences
    mt = mt.explode_rows(mt.sortedTranscriptConsequences)
    return mt.annotate_rows(
        info=mt.info.annotate(gene_id=mt.sortedTranscriptConsequences.gene_id)
    )


def annotate_sv1(
    mt: hl.MatrixTable, green_expression: hl.SetExpression
) -> hl.MatrixTable:
    """
    Annotate SVs with the SV1 category
    Rare, LOF, in a green gene

    Args:
        mt (hl.MatrixTable): the input MT
        green_expression (hl.SetExpression): the set of green genes

    Returns:
        hl.MatrixTable: the annotated MT
    """

    # filter again to LOF in green genes
    return mt.annotate_rows(
        info=mt.info.annotate(
            categorybooleansv1=hl.if_else(
                (mt.sortedTranscriptConsequences.major_consequence == 'LOF')
                & (green_expression.contains(mt.sortedTranscriptConsequences.gene_id)),
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def filter_matrix_by_ac(
    mt: hl.MatrixTable, ac_threshold: float | None = 0.05
) -> hl.MatrixTable:
    """
    Remove variants with AC in joint-call over threshold
    We don't need to worry about minimum cohort size
    due to the minimum group size of GATK-SV

    Args:
        mt (hl.MatrixTable):
        ac_threshold (float): remove variants more common than this in JointCall
    Returns:
        MT with all common-in-this-JC variants removed
    """

    return mt.filter_rows(
        (mt.info.MALE_AF[0] <= ac_threshold) & (mt.info.FEMALE_AF[0] <= ac_threshold)
    )


def rearrange_filters(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Rearrange the variantId MT, and remove filter-failing variants

    Args:
        mt ():
    """
    mt = mt.filter_rows(hl.is_missing(mt.filters) | (mt.filters.length() == 0))
    return mt.annotate_rows(
        info=mt.info.annotate(
            variantId=mt.variantId,
        )
    )


def main(
    mt_path: str,
    panelapp_path: str,
    pedigree: str,
    vcf_out: str,
):
    """
    Read MT, filter, and apply category annotation
    Export as a VCF

    Args:
        mt_path (str): where to find vcf output
        panelapp_path ():
        pedigree ():
        vcf_out (str): where to write VCF out
    """

    # read the parsed panelapp data
    get_logger().info(f'Reading PanelApp data from {panelapp_path!r}')
    panelapp = read_json_from_path(panelapp_path)['genes']  # type: ignore

    # pull green and new genes from the panelapp data
    # new is not currently incorporated in this analysis
    green_expression, _new_expression = green_and_new_from_panelapp(panelapp)

    # initiate Hail with defined driver spec.
    get_logger().info(f'Starting Hail with reference genome {genome_build()}')
    # un-comment this to play locally
    # hl.init(default_reference=genome_build())
    init_batch(driver_cores=1, driver_memory='lowmem')

    # if we already generated the annotated output, load instead
    if not to_path(mt_path.rstrip('/') + '/').exists():
        raise FileExistsError(f'Input MatrixTable doesn\'t exist: {mt_path}')

    # read in the input data (annotated)
    mt = hl.read_matrix_table(mt_path)

    # subset to currently considered samples
    mt = subselect_mt_to_pedigree(mt, pedigree=pedigree)

    # apply blanket filters
    mt = filter_matrix_by_ac(mt)
    mt = filter_matrix_by_af(mt)
    mt = rearrange_filters(mt)

    # pre-filter the MT and rearrange fields for export
    mt = restructure_mt_by_gene(mt)

    # label some SVs
    mt = annotate_sv1(mt, green_expression)
    # add further category annotations here

    # filter to labelled entries
    mt = mt.filter_rows(mt.info.categorybooleansv1 == ONE_INT)

    # now write that badboi
    hl.export_vcf(mt, vcf_out, tabix=True)


if __name__ == '__main__':

    # general CLI identical to the small variant version
    parser = ArgumentParser()
    parser.add_argument('--mt', required=True, help='path to input MT')
    parser.add_argument('--panelapp', type=str, required=True, help='panelapp JSON')
    parser.add_argument('--pedigree', type=str, required=True, help='Cohort Pedigree')
    parser.add_argument('--vcf_out', help='Where to write the VCF', required=True)
    args = parser.parse_args()

    logger = get_logger(__file__)
    logger.info('Running Hail filtering process for SVs')

    main(
        mt_path=args.mt,
        panelapp_path=args.panelapp,
        pedigree=args.pedigree,
        vcf_out=args.vcf_out,
    )
