"""
General purpose dumpster for hail methods common to
reanalysis & additional findings
"""

import logging

import hail as hl


def checkpoint_and_repartition(
    mt: hl.MatrixTable,
    checkpoint_root: str,
    checkpoint_num: int = 0,
    extra_logging: str | None = '',
) -> hl.MatrixTable:
    """
    uses an estimate of row size to inform the repartitioning of a MT
    aiming for a target partition size of ~10MB
    Kat's thread:
    https://discuss.hail.is/t/best-way-to-repartition-heavily-filtered-matrix-tables/2140
    :param mt:
    :param checkpoint_root:
    :param checkpoint_num:
    :param extra_logging: any additional context
    :return: repartitioned, post-checkpoint matrix
    """
    checkpoint_extended = f'{checkpoint_root}_{checkpoint_num}'
    logging.info(f'Checkpointing MT to {checkpoint_extended}')
    mt = mt.checkpoint(checkpoint_extended, overwrite=True)

    # estimate partitions; fall back to 1 if low row count
    current_rows = mt.count_rows()
    partitions = current_rows // 200000 or 1

    logging.info(
        f'Re-partitioning {current_rows} into {partitions} partitions {extra_logging}'
    )

    return mt.repartition(n_partitions=partitions, shuffle=True)


def filter_to_population_rare(mt: hl.MatrixTable, thresh: float) -> hl.MatrixTable:
    """
    run the rare filter, using Gnomad Exomes and Genomes
    :param mt:
    :param thresh: the float threshold we're using for gnomad
    :return:
    """
    # gnomad exomes and genomes below threshold or missing
    # if missing they were previously replaced with 0.0
    return mt.filter_rows(
        (mt.info.gnomad_ex_af < thresh) & (mt.info.gnomad_af < thresh)
    )


def filter_matrix_by_ac(
    mt: hl.MatrixTable, ac_threshold: float | None = 0.01
) -> hl.MatrixTable:
    """
    if called, this method will remove all variants in the joint call where the
    AlleleCount as a proportion is higher than the provided threshold
    :param mt:
    :param ac_threshold:
    :return: reduced MatrixTable
    """
    return mt.filter_rows((mt.info.AC <= 5) | (mt.info.AC / mt.info.AN < ac_threshold))


def filter_on_quality_flags(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    filter MT to rows with 0 quality filters
    note: in Hail, PASS is represented as an empty set
    :param mt:
    """
    return mt.filter_rows(mt.filters.length() == 0)


def filter_to_well_normalised(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    single alt per row, no missing Alt
    :param mt:
    """
    return mt.filter_rows((hl.len(mt.alleles) == 2) & (mt.alleles[1] != '*'))


def fields_audit(mt: hl.MatrixTable, fields: dict) -> bool:
    """
    checks that the required fields are all present before continuing
    """
    problems = []
    for field_group, group_types in fields.items():
        if field_group not in mt.row_value:
            problems.append(f'{field_group}:missing')
        else:
            for annotation, datatype in group_types:
                if annotation in mt[field_group]:
                    if not isinstance(mt[field_group][annotation], datatype):
                        problems.append(
                            f'{annotation}:'
                            f'{datatype}/'
                            f'{type(mt[field_group][annotation])}'
                        )
                else:
                    problems.append(f'{annotation}:missing')
    if problems:
        for problem in problems:
            logging.error(f'MT field: \t{problem}')
        return False
    return True


def vep_tx_audit(
    mt: hl.MatrixTable, vep_required: list[tuple[str, hl.Expression]]
) -> bool:
    """
    check that the required VEP annotations are present
    True if the 'audit' passes (all required fields present)
    """

    problems = []
    # now the content of the transcript_consequences
    if 'vep' not in mt.row_value:
        problems.append('VEP:missing')
    elif 'transcript_consequences' not in mt.vep:
        problems.append('transcript_consequences:missing')
    else:
        fields_and_types = dict(mt.vep.transcript_consequences[0].items())
        for field, field_type in vep_required:
            if field in fields_and_types:
                if not isinstance(fields_and_types[field], field_type):
                    problems.append(
                        f'{field}:{field_type}/{type(fields_and_types[field])}'
                    )
            else:
                problems.append(f'{field}:missing')

    if problems:
        logging.error('VEP field: \n'.join(problems))
        return False
    return True


def split_rows_by_gene_and_filter_to_green(
    mt: hl.MatrixTable, green_genes: hl.SetExpression
) -> hl.MatrixTable:
    """
    splits each GeneId onto a new row, then filters any
    rows not annotating a Green PanelApp gene

    geneIds is ENSG
    :param mt:
    :param green_genes:
    """

    # split each gene onto a separate row
    # transforms 'geneIds' field from set to string
    mt = mt.explode_rows(mt.geneIds)

    # filter rows without a green gene (removes empty geneIds)
    mt = mt.filter_rows(green_genes.contains(mt.geneIds))

    # limit the per-row transcript consequences to those relevant to the single
    # gene now present on each row
    return mt.annotate_rows(
        vep=mt.vep.annotate(
            transcript_consequences=mt.vep.transcript_consequences.filter(
                lambda x: (mt.geneIds == x.gene_id)
                & ((x.biotype == 'protein_coding') | (x.mane_select.contains('NM')))
            )
        )
    )
