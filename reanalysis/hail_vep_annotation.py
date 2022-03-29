"""
Read, filter, annotate, classify, and write Genetic data
- read VCF into MT
- read PanelApp data through GCP client
- hard-filter (FILTERS, AC)
- annotate
- extract generic fields
- remove all rows and consequences not relevant to GREEN genes
- consequence filter
- remove all rows with no interesting consequences
- extract vep data into CSQ string(s)
- annotate with categories 1, 2, 3, and 4
- remove all un-categorised variants
- write as VCF

This doesn't include applying inheritance pattern filters
Categories applied here are treated as unconfirmed
"""

from typing import Any, Dict, Optional
import json
import logging
import sys

import click
import hail as hl

from cloudpathlib import AnyPath


def filter_matrix_by_ac(
    matrix_data: hl.MatrixTable, config: Dict[str, Any]
) -> hl.MatrixTable:
    """

    :param matrix_data:
    :param config:
    :return: reduced MatrixTable
    """

    # count the samples in the VCF, and use to decide whether to implement
    # 'common within this joint call' as a filter
    # if we reach the sample threshold, filter on AC
    if matrix_data.count_cols() >= config['min_samples_to_ac_filter']:
        matrix_data = matrix_data.filter_rows(
            matrix_data.info.AC / matrix_data.info.AN < config['ac_threshold']
        )
    return matrix_data


def filter_matrix_by_variant_attributes(
    matrix_data: hl.MatrixTable, vqsr_run: Optional[bool] = True
) -> hl.MatrixTable:
    """
    filter MT to rows with normalised, high quality variants
    Note - when reading data into a MatrixTable, the Filters column is modified
    - split into a set of all filters
    - PASS is removed
    i.e. an empty set is equal to PASS in a VCF

    filter conditions applied are dependent on whether VQSR was run
    if VQSR - allow for empty filters, or VQSR with AS_FS=PASS
    if not - require the variant filters to be empty
    :param matrix_data:
    :param vqsr_run: if True, we
    :return:
    """
    vqsr_set = hl.literal({'VQSR'})
    pass_string = hl.literal('PASS')

    if vqsr_run:

        # hard filter for quality; assuming data is well normalised in pipeline
        matrix_data = matrix_data.filter_rows(
            (
                (matrix_data.filters.length() == 0)
                | (
                    (matrix_data.filters == vqsr_set)
                    & (matrix_data.info.AS_FilterStatus == pass_string)
                )
            )
        )

    # otherwise strictly enforce FILTERS==PASS, i.e. empty set
    else:
        matrix_data = matrix_data.filter_rows(matrix_data.filters.length() == 0)

    # normalised variants check
    # prior to annotation, variants in the MatrixTable representation are removed where:
    # - more than two alleles are present (ref and alt)
    #   - prior to annotation, the variant data must be decomposed to split all alt.
    #     alleles onto a separate row, with the corresponding sample genotypes
    # - alternate allele called is missing (*)
    matrix_data = matrix_data.filter_rows(
        (hl.len(matrix_data.alleles) == 2) & (matrix_data.alleles[1] != '*')
    )
    return matrix_data


def informed_repartition(
    matrix: hl.MatrixTable,
    temporary_path: str,
    post_annotation: Optional[bool] = False,
    small_fragments: Optional[bool] = False,
):
    """
    uses an estimate of row size to inform the repartitioning of a MT
    aiming for a target partition size of ~10MB
    post-annotation rows are estimated at ~5kB
        - a recursive sys.getsizeof-like guess suggested ~140 Bytes :/
        - writing a row to text was closer to 20kB :/
    pre-annotation rows are assumed to be substantially smaller
    throwing 400k rows into pre-annotation blocks, 150k into post

    Kat's thread:
    https://discuss.hail.is/t/best-way-to-repartition-heavily-filtered-matrix-tables/2140

    :param matrix:
    :param temporary_path:
    :param post_annotation:
    :param small_fragments:
    :return: repartitioned matrix
    """

    # calculate partitions, falling back to 1 partition if size is too small
    current_rows = matrix.count_rows()
    logging.info(f'Rows prior to repartition: {current_rows}')
    if post_annotation:
        partitions = current_rows // 200000 or 1
    elif small_fragments:
        partitions = current_rows // 800 or 1
    else:
        partitions = current_rows // 100000 or 1
    logging.info(f'Target number of partitions: {partitions}')

    # repartition with the specified # partitions
    matrix.repartition(n_partitions=partitions, shuffle=True)

    # a quick write to a temp path, and a read from the same
    return matrix.checkpoint(temporary_path, overwrite=True)


@click.command()
@click.option('--mt_inpath', help='path to the matrix table to ingest')
@click.option('--config_path', help='path to a config dict')
@click.option(
    '--mt_outpath',
    help='path to export annotated MT to',
    required=True,
)
def main(mt_inpath: str, config_path: str, mt_outpath: Optional[str] = None):
    """
    Read the MT from disk
    Do filtering and class annotation
    Export as a VCF

    trial:
    repartition estimate - each row is approx 7kb post annotation,
    to make a 10MB partition this is ~130k rows

    :param mt_inpath: path to the MT directory
    :param config_path: path to the config json
    :param mt_outpath:
    """

    # get the run configuration JSON
    logging.info(f'Reading config dict from "{config_path}"')
    with open(AnyPath(config_path), encoding='utf-8') as handle:
        config_dict = json.load(handle)

    # find the config area specific to hail operations
    hail_config = config_dict.get('filter')

    logging.info(
        f'Starting Hail with reference genome "{hail_config.get("ref_genome")}"'
    )
    # initiate Hail with the specified reference
    hl.init(default_reference=hail_config.get('ref_genome'), quiet=True)

    logging.info(f'Loading MT from "{mt_inpath}"')
    matrix = hl.read_matrix_table(mt_inpath)
    logging.debug(f'Loaded pre-annotation MT, size: {matrix.count_rows()}')

    # hard filter entries in the MT prior to annotation
    logging.info('Hard filtering variants')
    matrix = filter_matrix_by_ac(matrix_data=matrix, config=hail_config)
    matrix = filter_matrix_by_variant_attributes(matrix_data=matrix)
    logging.info(f'Post-Hard filtering size: {matrix.count_rows()}')
    # informed_repartition(matrix, post_annotation=False, temporary_path=mt_tmp)

    # re-annotate using VEP
    logging.info('Annotating variants')
    matrix = hl.vep(matrix, config='file:///vep_data/vep-gcloud.json')

    # if a path is provided, dump the MT
    matrix.write(mt_outpath, overwrite=True)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )
    main()  # pylint: disable=E1120
