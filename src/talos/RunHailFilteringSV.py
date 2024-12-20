"""
A Hail filtering process for labelling analysis-relevant SVs
Initially this will only contain a single category

CategoryBooleanSV1:
- rare
- predicted LoF in a listed gene
"""

from argparse import ArgumentParser

import hail as hl

from talos.config import config_retrieve
from talos.models import PanelApp
from talos.RunHailFiltering import MISSING_INT, ONE_INT, green_from_panelapp, subselect_mt_to_pedigree
from talos.static_values import get_logger
from talos.utils import read_json_from_path


def filter_matrix_by_af(mt: hl.MatrixTable, af_threshold: float = 0.03) -> hl.MatrixTable:
    """
    Filter a MatrixTable on AF, allow AF to be missing

    Args:
        mt (hl.MatrixTable): the input MT
        af_threshold (float): filtering threshold in gnomad v2.1

    Returns:
        same MT, with common variants removed
    """

    return mt.filter_rows(hl.or_else(mt.info['gnomad_v2.1_sv_AF'], MISSING_INT) < af_threshold)


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
    return mt.annotate_rows(info=mt.info.annotate(gene_id=mt.sortedTranscriptConsequences.gene_id))


def annotate_sv1(mt: hl.MatrixTable, green_expression: hl.SetExpression) -> hl.MatrixTable:
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
            ),
        ),
    )


def filter_matrix_by_ac(mt: hl.MatrixTable, ac_threshold: float | None = 0.03) -> hl.MatrixTable:
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

    return mt.filter_rows((mt.info.MALE_AF[0] <= ac_threshold) & (mt.info.FEMALE_AF[0] <= ac_threshold))


def rearrange_variant_id(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Rearrange the variantId MT
    For now we do not filter out failing variant sites, we'll reserve this for per-sample quality failures

    Args:
        mt ():
    """
    return mt.annotate_rows(info=mt.info.annotate(variantId=mt.variantId))


def fix_hemi_calls(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Hail's MT -> VCF export doesn't handle hemizygous calls
    adjust the relevant single allele calls to a biallelic representation
    going with Hom-Alt/Hom-Ref

    if GT == 1, recast as [1, 1]
    if GT == 0, recast as [0, 0]

    Args:
        mt ():
    """

    return mt.annotate_entries(
        GT=hl.if_else(
            mt.GT.is_diploid(),
            mt.GT,
            hl.if_else(mt.GT.is_non_ref(), hl.call(1, 1), hl.call(0, 0)),
        ),
    )


def cli_main():
    """
    main method wrapper for console script execution
    """
    parser = ArgumentParser()
    parser.add_argument('--input', required=True, help='path to input MT')
    parser.add_argument('--panelapp', type=str, required=True, help='GeneratePanelData JSON')
    parser.add_argument('--pedigree', type=str, required=True, help='Cohort Pedigree')
    parser.add_argument('--output', help='Where to write the VCF', required=True)
    args = parser.parse_args()
    main(mt_path=args.input, panelapp_path=args.panelapp, pedigree=args.pedigree, vcf_out=args.output)


def main(mt_path: str, panelapp_path: str, pedigree: str, vcf_out: str):
    """
    Read MT, filter, and apply category annotation
    Export as a VCF

    Args:
        mt_path (str): where to find vcf output
        panelapp_path ():
        pedigree ():
        vcf_out (str): where to write VCF out
    """
    get_logger(__file__).info(
        r"""Welcome To
███████████   █████████   █████          ███████     █████████
█   ███   █  ███     ███   ███         ███     ███  ███     ███
    ███      ███     ███   ███        ███       ███ ███
    ███      ███████████   ███        ███       ███  █████████
    ███      ███     ███   ███        ███       ███         ███
    ███      ███     ███   ███      █  ███     ███  ███     ███
   █████    █████   █████ ███████████    ███████     █████████
        (SV style)""",
    )

    # read the parsed panelapp data
    get_logger().info(f'Reading PanelApp data from {panelapp_path!r}')
    panelapp = read_json_from_path(panelapp_path, return_model=PanelApp)
    if not isinstance(panelapp, PanelApp):
        raise TypeError(f'PanelApp was not a PanelApp object: {panelapp}')

    # pull green and new genes from the panelapp data
    # new is not currently incorporated in this analysis
    green_expression = green_from_panelapp(panelapp)

    # initiate Hail in local cluster mode
    number_of_cores = config_retrieve(['RunHailFiltering', 'cores', 'sv'], 2)
    get_logger().info(f'Starting Hail with reference genome GRCh38, as a {number_of_cores} core local cluster')
    hl.context.init_spark(master=f'local[{number_of_cores}]', quiet=True)
    hl.default_reference('GRCh38')

    # read in the input data (annotated) Not checking for existence; if it fails, it fails
    mt = hl.read_matrix_table(mt_path)

    # subset to currently considered samples
    mt = subselect_mt_to_pedigree(mt, pedigree=pedigree)

    # apply blanket filters
    mt = filter_matrix_by_ac(mt, ac_threshold=config_retrieve(['RunHailFiltering', 'callset_af_sv_recessive']))
    mt = filter_matrix_by_af(mt, af_threshold=config_retrieve(['RunHailFiltering', 'callset_af_sv_recessive']))
    mt = rearrange_variant_id(mt)

    # pre-filter the MT and rearrange fields for export
    mt = restructure_mt_by_gene(mt)

    # label some SVs
    mt = annotate_sv1(mt, green_expression)

    # add further category annotations here

    # filter to labelled entries
    mt = mt.filter_rows(mt.info.categorybooleansv1 == ONE_INT)

    # Hail's MT -> VCF export doesn't handle hemizygous calls
    mt = fix_hemi_calls(mt)

    # now write that badboi
    hl.export_vcf(mt, vcf_out, tabix=True)


if __name__ == '__main__':
    cli_main()
