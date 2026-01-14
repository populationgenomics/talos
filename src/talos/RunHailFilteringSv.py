"""
A Hail filtering process for labelling analysis-relevant SVs
Initially this will only contain a single category
This expects data annotated by GATK's SVAnnotate tool

CategoryBooleanSV1:
- rare
- predicted LoF in a listed gene
"""

from argparse import ArgumentParser

from loguru import logger
from mendelbrot.pedigree_parser import PedigreeParser

import hail as hl

from talos.config import config_retrieve
from talos.models import PanelApp
from talos.RunHailFiltering import MISSING_INT, ONE_INT, green_from_panelapp, subselect_mt_to_pedigree
from talos.utils import read_json_from_path

GNOMAD_POP = config_retrieve(['RunHailFilteringSv', 'gnomad_population'], 'gnomad_v4.1')


def read_and_filter_mane_json(mane_json: str) -> hl.dict:
    """
    Read the MANE JSON and filter it to the relevant fields
    Args:
        mane_json ():

    Returns:

    """

    json_dict = read_json_from_path(mane_json)

    return hl.literal({entry['symbol']: entry['ensg'] for entry in json_dict.values()})


def rearrange_annotations(mt: hl.MatrixTable, gene_mapping: hl.dict) -> hl.MatrixTable:
    """
    Rearrange the annotations in the MT to be more easily accessible
    Args:
        mt ():
        gene_mapping (): gene symbol to gene ID mapping

    Returns:
        same MT, with shifted annotations
    """

    # accept, but don't force, this GATK-SV field
    if 'ALGORITHMS' not in mt.info:
        mt = mt.annotate_rows(
            info=mt.info.annotate(
                ALGORITHMS=['gCNV'],
            ),
        )

    # embed these attributes if missing, but give null values
    for sv_attribute in ('STATUS', 'CHR2', 'END2'):
        if sv_attribute not in mt.info:
            mt = mt.annotate_rows(
                info=mt.info.annotate(
                    **{sv_attribute: hl.missing(hl.tstr)},
                ),
            )

    # trying to sustain backwards compatibility with a changing GATK-SV/gCNV pipeline combination
    if 'AF_MALE' in mt.info:
        mt = mt.annotate_rows(
            info=mt.info.annotate(
                male_af=mt.info.AF_MALE,
                female_af=mt.info.AF_FEMALE,
            ),
        )
    else:
        mt = mt.annotate_rows(
            info=mt.info.annotate(
                male_af=mt.info.MALE_AF,
                female_af=mt.info.FEMALE_AF,
            ),
        )

    mt = mt.annotate_rows(
        info=hl.struct(
            AC=mt.info.AC,
            AF=mt.info.AF,
            AN=mt.info.AN,
            algorithms=mt.info.ALGORITHMS,
            gnomad_sv_ID=mt.info[f'{GNOMAD_POP}_sv_SVID'],
            gnomad_sv_AF=mt.info[f'{GNOMAD_POP}_sv_AF'],
            lof=hl.set(mt.info['PREDICTED_LOF']),
            n_het=mt.info.N_HET,
            n_homalt=mt.info.N_HOMALT,
            svlen=mt.info.SVLEN,
            svtype=mt.info.SVTYPE,
            status=mt.info.STATUS,
            end=mt.info.END,
            chr2=mt.info.CHR2,
            end2=mt.info.END2,
            male_af=mt.info.male_af,
            female_af=mt.info.female_af,
        ),
    )

    # match the symbols to gene IDs
    return mt.annotate_rows(
        info=mt.info.annotate(
            lof_ensg=hl.set(hl.map(lambda gene: gene_mapping.get(gene, gene), mt.info.lof)),
            # this is so we can explode it out later, whilst keeping the full list
            gene_id=hl.set(hl.map(lambda gene: gene_mapping.get(gene, gene), mt.info.lof)),
        ),
    )


def filter_matrix_by_af(mt: hl.MatrixTable, af_threshold: float = 0.03) -> hl.MatrixTable:
    """
    Filter a MatrixTable on AF, allow AF to be missing

    Args:
        mt (hl.MatrixTable): the input MT
        af_threshold (float): filtering threshold in gnomad v4.1

    Returns:
        same MT, with common variants removed
    """

    return mt.filter_rows(
        hl.or_else(
            mt.info.gnomad_sv_AF,
            MISSING_INT,
        )
        < af_threshold,
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

    return mt.filter_rows(
        (mt.info.male_af[0] <= ac_threshold) & (mt.info.female_af[0] <= ac_threshold),
    )


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
    parser.add_argument(
        '--input',
        required=True,
        help='path to input VCF',
    )
    parser.add_argument(
        '--panelapp',
        required=True,
        help='GeneratePanelData JSON',
    )
    parser.add_argument(
        '--mane_json',
        required=True,
        help='JSON for gene~symbol mapping',
    )
    parser.add_argument(
        '--pedigree',
        required=True,
        help='Cohort Pedigree',
    )
    parser.add_argument(
        '--output',
        help='Where to write the VCF',
        required=True,
    )
    args = parser.parse_args()
    main(
        vcf_path=args.input,
        panelapp_path=args.panelapp,
        mane_json=args.mane_json,
        pedigree=args.pedigree,
        vcf_out=args.output,
    )


def main(vcf_path: str, panelapp_path: str, mane_json: str, pedigree: str, vcf_out: str):
    """
    Read MT, filter, and apply category annotation
    Export as a VCF

    Args:
        vcf_path (str): where to find vcf output
        panelapp_path ():
        mane_json ():
        pedigree ():
        vcf_out (str): where to write VCF out
    """
    logger.info(
        r"""Welcome To
 ███████████   █████████   █████          ███████     █████████
 █   ███   █  ███     ███   ███         ███     ███  ███     ███
     ███      ███     ███   ███        ███       ███ ███
     ███      ███████████   ███        ███       ███  █████████
     ███      ███     ███   ███        ███       ███         ███
     ███      ███     ███   ███      █  ███     ███  ███     ███
    █████    █████   █████ ███████████    ███████     █████████
        """,
    )

    # get the Gene-Symbol mapping dict
    gene_id_mapping = read_and_filter_mane_json(mane_json)

    # read the parsed panelapp data
    logger.info(f'Reading PanelApp data from {panelapp_path!r}')
    panelapp = read_json_from_path(panelapp_path, return_model=PanelApp)
    if not isinstance(panelapp, PanelApp):
        raise TypeError(f'PanelApp was not a PanelApp object: {panelapp}')

    # pull green and new genes from the panelapp data
    # new is not currently incorporated in this analysis
    green_expression = green_from_panelapp(panelapp)

    # initiate Hail in local cluster mode
    logger.info('Starting Hail with reference genome GRCh38, as a local cluster')
    hl.context.init_spark(master='local[*]')
    hl.default_reference('GRCh38')

    # read the VCF in as a MatrixTable, and checkpoint it locally
    mt = hl.import_vcf(
        vcf_path,
        reference_genome='GRCh38',
        skip_invalid_loci=True,
        force_bgz=True,
    ).checkpoint(output='temporary.mt')

    # parse the pedigree into an object
    pedigree_data = PedigreeParser(pedigree)

    # reduce cohort to singletons, if the config says so
    if config_retrieve('singletons', False):
        logger.info('Reducing pedigree to affected singletons only')
        pedigree_data.set_participants(pedigree_data.as_singletons())
        pedigree_data.set_participants(pedigree_data.get_affected_members())

    # subset to currently considered samples
    mt = subselect_mt_to_pedigree(mt, ped_samples=pedigree_data.get_all_sample_ids())

    # remove filtered variants
    mt = mt.filter_rows(hl.is_missing(mt.filters) | (mt.filters.length() == 0))

    logger.info(f'Loaded {mt.count_rows()} rows and {mt.count_cols()} columns, in {mt.n_partitions()} partitions')

    # drop rows with no LOF consequences
    mt = mt.filter_rows(hl.len(mt.info.PREDICTED_LOF) > 0)

    # rearrange the annotations
    mt = rearrange_annotations(mt, gene_id_mapping)

    # apply blanket filters
    ac_threshold = config_retrieve(['RunHailFiltering', 'callset_af_sv_recessive'])
    mt = filter_matrix_by_ac(mt, ac_threshold=ac_threshold)
    mt = filter_matrix_by_af(mt, af_threshold=ac_threshold)

    mt = mt.explode_rows(mt.info.gene_id)

    # hard filter remaining rows to PanelApp green genes
    mt = mt.filter_rows(green_expression.contains(mt.info.gene_id))

    # everything left is `SV1`
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            categorybooleansv1=ONE_INT,
        ),
    )

    # Hail's MT -> VCF export doesn't handle hemizygous calls
    mt = fix_hemi_calls(mt)

    # now write that badboi
    hl.export_vcf(
        mt,
        vcf_out,
        tabix=True,
    )


if __name__ == '__main__':
    cli_main()
