"""
Read, filter, annotate, classify, and write Genetic data
- read MT
- read PanelApp data through GCP client
- hard-filter (FILTERS, AC)
- extract generic fields
- remove all rows and consequences not relevant to GREEN genes
- extract vep data into CSQ string(s)
- annotate with categories 1, 2, 3, 4, 5, and Support
- remove un-categorised variants
- write as VCF
This doesn't include applying inheritance pattern filters
Categories applied here are treated as unconfirmed
"""

from typing import Any
import json
import logging
import sys
from argparse import ArgumentParser

import hail as hl

from cloudpathlib import AnyPath
from cpg_utils.hail_batch import init_batch


# set some Hail constants
MISSING_INT = hl.int32(0)
MISSING_FLOAT_LO = hl.float64(0.0)
MISSING_FLOAT_HI = hl.float64(1.0)
MISSING_STRING = hl.str('missing')
ONE_INT = hl.int32(1)
BENIGN = hl.str('benign')
CONFLICTING = hl.str('conflicting')
LOFTEE_HC = hl.str('HC')
PATHOGENIC = hl.str('pathogenic')


def filter_matrix_by_ac(
    matrix: hl.MatrixTable, ac_threshold: float | None = 0.1
) -> hl.MatrixTable:
    """
    if called, this method will remove all variants in the joint call where the
    AlleleCount as a proportion is higher than the provided threshold
    :param matrix:
    :param ac_threshold:
    :return: reduced MatrixTable
    """
    return matrix.filter_rows(matrix.info.AC / matrix.info.AN < ac_threshold)


def filter_on_quality_flags(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    filter MT to rows with 0 quality filters
    note: in Hail, PASS is represented as an empty set
    :param matrix:
    """
    return matrix.filter_rows(matrix.filters.length() == 0)


def filter_to_well_normalised(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    single alt per row, no missing Alt
    :param matrix:
    """
    return matrix.filter_rows(
        (hl.len(matrix.alleles) == 2) & (matrix.alleles[1] != '*')
    )


def filter_by_ab_ratio(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    filters HomRef, Het, and HomAlt by appropriate AB ratio bins
    :param matrix:
    """
    ab = matrix.AD[1] / hl.sum(matrix.AD)
    return matrix.filter_entries(
        (matrix.GT.is_hom_ref() & (ab <= 0.15))
        | (matrix.GT.is_het() & (ab >= 0.25) & (ab <= 0.75))
        | (matrix.GT.is_hom_var() & (ab >= 0.85))
    )


def annotate_category_1(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    applies the Category1 annotation (1 or 0) as appropriate
    semi-rare in Gnomad
    at least one Clinvar star
    contains pathogenic and not conflicting; doesn't contain benign
    :param matrix:
    :return: same Matrix, with additional field per variant
    """

    return matrix.annotate_rows(
        info=matrix.info.annotate(
            Category1=hl.if_else(
                (matrix.info.clinvar_stars > 0)
                & (matrix.info.clinvar_sig.lower().contains(PATHOGENIC))
                & ~(matrix.info.clinvar_sig.lower().contains(CONFLICTING))
                & ~(matrix.info.clinvar_sig.lower().contains(BENIGN)),
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def annotate_category_2(
    matrix: hl.MatrixTable, config: dict[str, Any], new_genes: hl.SetExpression
) -> hl.MatrixTable:
    """
    - Gene is new in PanelApp
    - Clinvar contains pathogenic, or
    - Critical protein consequence on at least one transcript
    - High in silico consequence
    :param matrix:
    :param config:
    :param new_genes: the new genes in this panelapp content
    :return: same Matrix, with additional field per variant
    """

    critical_consequences = hl.set(config.get('critical_csq'))

    # check for new - if new, allow for in silico, CSQ, or clinvar to confirm
    return matrix.annotate_rows(
        info=matrix.info.annotate(
            Category2=hl.if_else(
                (new_genes.contains(matrix.geneIds))
                & (
                    (
                        hl.len(
                            matrix.vep.transcript_consequences.filter(
                                lambda x: hl.len(
                                    critical_consequences.intersection(
                                        hl.set(x.consequence_terms)
                                    )
                                )
                                > 0
                            )
                        )
                        > 0
                    )
                    | (
                        (matrix.info.clinvar_sig.lower().contains(PATHOGENIC))
                        & ~(matrix.info.clinvar_sig.lower().contains(CONFLICTING))
                        & ~(matrix.info.clinvar_sig.lower().contains(BENIGN))
                    )
                    | (
                        (matrix.info.cadd > config['in_silico']['cadd'])
                        | (matrix.info.revel > config['in_silico']['revel'])
                    )
                ),
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def annotate_category_3(
    matrix: hl.MatrixTable, config: dict[str, Any]
) -> hl.MatrixTable:
    """
    applies the Category3 flag where appropriate
    - Critical protein consequence on at least one transcript
    - either predicted NMD or
    - any star Pathogenic or Likely_pathogenic in Clinvar
    :param matrix:
    :param config:
    :return:
    """

    critical_consequences = hl.set(config.get('critical_csq'))

    # First check if we have any HIGH consequences
    # then explicitly link the LOFTEE check with HIGH consequences
    # OR allow for a pathogenic ClinVar, any Stars
    return matrix.annotate_rows(
        info=matrix.info.annotate(
            Category3=hl.if_else(
                (
                    hl.len(
                        matrix.vep.transcript_consequences.filter(
                            lambda x: (
                                hl.len(
                                    critical_consequences.intersection(
                                        hl.set(x.consequence_terms)
                                    )
                                )
                                > 0
                            )
                        )
                    )
                    > 0
                )
                & (
                    (
                        hl.len(
                            matrix.vep.transcript_consequences.filter(
                                lambda x: (
                                    hl.len(
                                        critical_consequences.intersection(
                                            hl.set(x.consequence_terms)
                                        )
                                    )
                                    > 0
                                )
                                & ((x.lof == LOFTEE_HC) | (hl.is_missing(x.lof)))
                            )
                        )
                        > 0
                    )
                    | (
                        (matrix.info.clinvar_sig.lower().contains(PATHOGENIC))
                        & ~(matrix.info.clinvar_sig.lower().contains(CONFLICTING))
                        & ~(matrix.info.clinvar_sig.lower().contains(BENIGN))
                    )
                ),
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def filter_by_consequence(
    matrix: hl.MatrixTable, config: dict[str, Any]
) -> hl.MatrixTable:
    """
    - reduce the per-row transcript consequences to a limited group
    - reduce the rows to ones where there are remaining tx consequences
    :param matrix:
    :param config: dictionary content relating to hail
    :return: reduced matrix
    """

    # at time of writing this is VEP HIGH + missense_variant
    # update without updating the dictionary content
    high_csq = set(config.get('critical_csq', [])).union(
        set(config.get('additional_consequences', []))
    )

    # overwrite the consequences with an intersection against a limited list
    matrix = matrix.annotate_rows(
        vep=matrix.vep.annotate(
            transcript_consequences=matrix.vep.transcript_consequences.filter(
                lambda x: hl.len(hl.set(x.consequence_terms).intersection(high_csq)) > 0
            )
        )
    )

    # filter out rows with no tx consequences left, and no splice cat. assignment
    return matrix.filter_rows(
        (hl.len(matrix.vep.transcript_consequences) > 0) & (matrix.info.Category5 == 0)
    )


def annotate_category_4(
    matrix: hl.MatrixTable, config: dict[str, Any], plink_family_file: str
) -> hl.MatrixTable:
    """
    Category based on de novo MOI, restricted to a group of consequences
    Initial implementation didn't restrict by tx consequence, which meant we
    found 10s of variants per sample, all intronic variants
    This category focuses de novo on MOI within family structures
    Instead, we are leveraging the inbuilt hl.de_novo functionality to:
    - identify all likely de novo inherited variation
    - collect samples for each variant (expect n=1 per de novo call, but not guaranteed)
    - compress list of proband IDs as a String
    - apply string back to the MatrixTable, else 'missing' if missing
    :param matrix:
    :param config:
    :param plink_family_file: path to a pedigree in PLINK format
    :return:
    """

    logging.info('running de novo search')

    de_novo_matrix = filter_by_consequence(matrix, config)

    # read pedigree from the specified file
    pedigree = hl.Pedigree.read(plink_family_file)

    # identify likely de novo calls
    # use the AFs already embedded in the MT as the prior population frequencies
    dn_table = hl.de_novo(
        de_novo_matrix, pedigree, pop_frequency_prior=de_novo_matrix.info.gnomad_af
    )
    dn_table = dn_table.filter(dn_table.confidence == 'HIGH')

    # re-key the table by locus,alleles, removing the sampleID from the compound key
    dn_table = dn_table.key_by(dn_table.locus, dn_table.alleles)

    # we only require the key (locus, alleles) and the sample ID
    # select to remove other fields, then collect per-key into Array of Structs
    dn_table = dn_table.select(dn_table.id).collect_by_key()

    # collect all sample IDs per locus, and squash into a String Array
    # delimit to compress that Array into single Strings
    dn_table = dn_table.annotate(
        values=hl.delimit(hl.map(lambda x: x.id, dn_table.values), ',')
    )

    # log the number of variants found this way
    logging.info(f'{dn_table.count()} variants showed de novo inheritance')

    # annotate those values as a flag if relevant, else 'missing'
    return matrix.annotate_rows(
        info=matrix.info.annotate(
            Category4=hl.or_else(dn_table[matrix.row_key].values, MISSING_STRING)
        )
    )


def annotate_category_5(
    matrix: hl.MatrixTable, config: dict[str, Any]
) -> hl.MatrixTable:
    """
    :param matrix:
    :param config:
    :return: same Matrix, with additional field per variant
    """

    return matrix.annotate_rows(
        info=matrix.info.annotate(
            Category5=hl.if_else(
                matrix.info.splice_ai_delta >= config['spliceai_full'],
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def annotate_category_support(
    matrix: hl.MatrixTable, config: dict[str, Any]
) -> hl.MatrixTable:
    """
    Background class based on in silico annotations
    - rare in Gnomad, and
    - CADD & REVEL above threshold (switched to consensus), or
    - Massive cross-tool consensus
    - polyphen and sift are evaluated per-consequence
    The intention is that this class will never be the sole reason to treat a variant
    as interesting, but can be used as a broader class with a lower barrier to entry
    Goal is to have a broader set of variants to find Comp-Hets.
    This support category can be expanded to encompass other scenarios that qualify
    variants as 'of second-hit interest only'
    :param matrix:
    :param config:
    :return:
    """

    return matrix.annotate_rows(
        info=matrix.info.annotate(
            CategorySupport=hl.if_else(
                (
                    (matrix.info.cadd > config['in_silico'].get('cadd'))
                    & (matrix.info.revel > config['in_silico'].get('revel'))
                )
                | (
                    (
                        matrix.vep.transcript_consequences.any(
                            lambda x: hl.or_else(x.sift_score, MISSING_FLOAT_HI)
                            <= config['in_silico'].get('sift')
                        )
                    )
                    & (
                        matrix.vep.transcript_consequences.any(
                            lambda x: hl.or_else(x.polyphen_score, MISSING_FLOAT_LO)
                            >= config['in_silico'].get('polyphen')
                        )
                    )
                    & (
                        (matrix.info.mutationtaster.contains('D'))
                        | (matrix.info.mutationtaster == 'missing')
                    )
                ),
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def transform_variant_string(locus_details: hl.Struct) -> str:
    """
    takes an object
    Struct(
        locus=Locus(
            contig='chr1',
            position=10,
            reference_genome='GRCh38'
        ),
        alleles=['GC', 'G'],
        category_4_only=0
    )
    transform into simplified 1-10-GC-G
    drop the category_4_only attribute
    :param locus_details:
    :return:
    """
    return '-'.join(
        [
            locus_details.locus.contig.replace('chr', ''),
            str(locus_details.locus.position),
            *locus_details.alleles,
        ]
    )


def filter_to_population_rare(
    matrix: hl.MatrixTable, config: dict[str, Any]
) -> hl.MatrixTable:
    """
    run the rare filter, using Gnomad & exac
    :param matrix:
    :param config:
    :return:
    """
    # exac and gnomad must be below threshold or missing
    # if missing they were previously replaced with 0.0
    # could also extend this filter to include max gnomad Homs
    return matrix.filter_rows(
        (matrix.info.exac_af < config['af_semi_rare'])
        & (matrix.info.gnomad_af < config['af_semi_rare'])
    )


def split_rows_by_gene_and_filter_to_green(
    matrix: hl.MatrixTable, green_genes: hl.SetExpression
) -> hl.MatrixTable:
    """
    splits each GeneId onto a new row, then filters any
    rows not annotating a Green PanelApp gene
    :param matrix:
    :param green_genes:
    """

    # split each gene onto a separate row
    # transforms 'geneIds' field from set to string
    matrix = matrix.explode_rows(matrix.geneIds)

    # filter rows without a green gene (removes empty geneIds)
    matrix = matrix.filter_rows(green_genes.contains(matrix.geneIds))

    # limit the per-row transcript consequences to those relevant to the single
    # gene now present on each row
    matrix = matrix.annotate_rows(
        vep=matrix.vep.annotate(
            transcript_consequences=matrix.vep.transcript_consequences.filter(
                lambda x: (matrix.geneIds == x.gene_id)
                & ((x.biotype == 'protein_coding') | (x.mane_select.contains('NM')))
            )
        )
    )

    return matrix


def vep_struct_to_csq(
    vep_expr: hl.expr.StructExpression, csq_fields: list[str]
) -> hl.expr.ArrayExpression:
    """
    Taken shamelessly from the gnomad library source code
    Given a VEP Struct, returns an array of VEP VCF CSQ strings
    (1 per csq in the struct).
    Fields & order correspond to those in `csq_fields`, corresponding to the
    VCF header that is required to interpret the VCF CSQ INFO field.
    Order is flexible & all fields in the default value are supported.
    These fields are formatted in the same way that their VEP CSQ counterparts are.
    :param vep_expr: The input VEP Struct
    :param csq_fields: ordered list of lower-case fields to include in the CSQ
    :return: The corresponding CSQ string(s)
    """

    def get_csq_from_struct(
        element: hl.expr.StructExpression, feature_type: str
    ) -> hl.expr.StringExpression:

        # Most fields are 1-1, just lowercase
        fields = dict(element)

        # Add general exceptions
        fields.update(
            {
                'allele': element.variant_allele,
                'consequence': hl.delimit(element.consequence_terms, delimiter='&'),
                'feature_type': feature_type,
                'feature': (
                    element.transcript_id
                    if 'transcript_id' in element
                    else element.regulatory_feature_id
                    if 'regulatory_feature_id' in element
                    else element.motif_feature_id
                    if 'motif_feature_id' in element
                    else ''
                ),
                'variant_class': vep_expr.variant_class,
                'canonical': hl.if_else(element.canonical == 1, 'YES', ''),
                'ensp': element.protein_id,
                'gene': element.gene_id,
                'symbol': element.gene_symbol,
                'symbol_source': element.gene_symbol_source,
                'cdna_position': hl.str(element.cdna_start)
                + hl.if_else(
                    element.cdna_start == element.cdna_end,
                    '',
                    '-' + hl.str(element.cdna_end),
                ),
                'cds_position': hl.str(element.cds_start)
                + hl.if_else(
                    element.cds_start == element.cds_end,
                    '',
                    '-' + hl.str(element.cds_end),
                ),
                'protein_position': hl.str(element.protein_start)
                + hl.if_else(
                    element.protein_start == element.protein_end,
                    '',
                    '-' + hl.str(element.protein_end),
                ),
                'sift': element.sift_prediction
                + '('
                + hl.format('%.3f', element.sift_score)
                + ')',
                'polyphen': element.polyphen_prediction
                + '('
                + hl.format('%.3f', element.polyphen_score)
                + ')',
                'mane_select': element.mane_select,
            }
        )

        return hl.delimit(
            [hl.or_else(hl.str(fields.get(f, '')), '') for f in csq_fields], '|'
        )

    csq = hl.empty_array(hl.tstr)
    csq = csq.extend(
        hl.or_else(
            vep_expr['transcript_consequences'].map(
                lambda x: get_csq_from_struct(x, feature_type='Transcript')
            ),
            hl.empty_array(hl.tstr),
        )
    )

    # prior filtering on consequence will make this caution unnecessary
    return hl.or_missing(hl.len(csq) > 0, csq)


def extract_annotations(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    pull out select fields which aren't per-consequence
    store these in INFO (required to be included in VCF export)
    replace with placeholder (least consequential) if empty
    e.g. most tools score 0, but for Sift 1 is least important
    :param matrix:
    :return: input matrix with annotations pulled into INFO
    """

    return matrix.annotate_rows(
        info=matrix.info.annotate(
            exac_af=hl.or_else(matrix.exac.AF, MISSING_FLOAT_LO),
            exac_ac_het=hl.or_else(matrix.exac.AC_Het, MISSING_INT),
            exac_ac_hom=hl.or_else(matrix.exac.AC_Hom, MISSING_INT),
            exac_ac_hemi=hl.or_else(matrix.exac.AC_Hemi, MISSING_INT),
            gnomad_ex_cov=hl.or_else(matrix.gnomad_exome_coverage, MISSING_FLOAT_LO),
            gnomad_ex_af=hl.or_else(matrix.gnomad_exomes.AF, MISSING_FLOAT_LO),
            gnomad_ex_an=hl.or_else(matrix.gnomad_exomes.AN, MISSING_INT),
            gnomad_ex_ac=hl.or_else(matrix.gnomad_exomes.AC, MISSING_INT),
            gnomad_ex_hom=hl.or_else(matrix.gnomad_exomes.Hom, MISSING_INT),
            gnomad_cov=hl.or_else(matrix.gnomad_genome_coverage, MISSING_FLOAT_LO),
            gnomad_af=hl.or_else(matrix.gnomad_genomes.AF, MISSING_FLOAT_LO),
            gnomad_an=hl.or_else(matrix.gnomad_genomes.AN, MISSING_INT),
            gnomad_ac=hl.or_else(matrix.gnomad_genomes.AC, MISSING_INT),
            gnomad_hom=hl.or_else(matrix.gnomad_genomes.Hom, MISSING_INT),
            splice_ai_delta=hl.or_else(matrix.splice_ai.delta_score, MISSING_FLOAT_LO),
            splice_ai_csq=hl.or_else(
                matrix.splice_ai.splice_consequence, MISSING_STRING
            ).replace(' ', '_'),
            revel=hl.float64(hl.or_else(matrix.dbnsfp.REVEL_score, '0.0')),
            cadd=hl.or_else(matrix.cadd.PHRED, MISSING_FLOAT_LO),
            clinvar_sig=hl.or_else(
                matrix.clinvar.clinical_significance, MISSING_STRING
            ),
            clinvar_stars=hl.or_else(matrix.clinvar.gold_stars, MISSING_INT),
            # these next 3 are per-transcript, with ';' to delimit
            # pulling these annotations into INFO with ';' to separate
            # will break INFO parsing for most tools
            mutationtaster=hl.or_else(
                matrix.dbnsfp.MutationTaster_pred, MISSING_STRING
            ).replace(';', ','),
            fathmm=hl.or_else(matrix.dbnsfp.FATHMM_pred, MISSING_STRING).replace(
                ';', ','
            ),
            metasvm=hl.or_else(matrix.dbnsfp.MetaSVM_pred, MISSING_STRING).replace(
                ';', ','
            ),
            phast_cons=hl.float64(
                hl.or_else(matrix.dbnsfp.phastCons100way_vertebrate, '0.0')
            ),
            gerp_rs=hl.float64(hl.or_else(matrix.dbnsfp.GERP_RS, '0.0')),
            eigen_phred=hl.or_else(matrix.eigen.Eigen_phred, MISSING_FLOAT_LO),
        )
    )


def filter_to_categorised(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    Filter to rows tagged with a class
    :param matrix:
    :return: input matrix, minus rows without Categories applied
    """
    return matrix.filter_rows(
        (matrix.info.Category1 == 1)
        | (matrix.info.Category2 == 1)
        | (matrix.info.Category3 == 1)
        | (matrix.info.Category4 != 'missing')
        | (matrix.info.Category5 == 1)
        | (matrix.info.CategorySupport == 1)
    )


def write_matrix_to_vcf(
    matrix: hl.MatrixTable, output_path: str, additional_header: str | None = None
):
    """
    write the remaining MatrixTable content to file as a VCF
    :param matrix:
    :param output_path: where to write
    :param additional_header: file containing any other lines to add into header
    """
    hl.export_vcf(
        matrix,
        output_path,
        append_to_header=additional_header,
        tabix=True,
    )


def green_and_new_from_panelapp(
    panel_data: dict[str, dict[str, str]]
) -> tuple[hl.SetExpression, hl.SetExpression]:
    """
    Pull all ENSGs from PanelApp data relating to Green Genes
    Also identify the subset of those genes which relate to NEW in panel
    :param panel_data:
    :return: two set expressions, green genes and new genes
    """

    # take all the green genes, remove the metadata
    green_genes = set(panel_data.keys()) - {'panel_metadata'}
    logging.info(f'Extracted {len(green_genes)} green genes')
    green_gene_set_expression = hl.literal(green_genes)

    new_genes = {gene for gene in green_genes if panel_data[gene].get('new')}
    logging.info(f'Extracted {len(new_genes)} NEW genes')
    new_gene_set_expression = hl.literal(new_genes)

    return green_gene_set_expression, new_gene_set_expression


def informed_repartition(matrix: hl.MatrixTable):
    """
    uses an estimate of row size to inform the repartitioning of a MT
    aiming for a target partition size of ~10MB
    Kat's thread:
    https://discuss.hail.is/t/best-way-to-repartition-heavily-filtered-matrix-tables/2140
    :param matrix:
    :return: repartitioned matrix
    """

    # estimate partitions; fall back to 1 if low row count
    current_rows = matrix.count_rows()
    partitions = current_rows // 200000 or 1

    logging.info(f'Re-partitioning {current_rows} into {partitions} partitions')
    return matrix.repartition(n_partitions=partitions, shuffle=True)


def main(
    mt_input: str,
    panelapp_path: str,
    config_path: str,
    plink_file: str,
    out_vcf: str,
    tmp_path: str,
):
    """
    Read the MT from disk
    Do filtering and class annotation
    Export as a VCF
    :param mt_input: path to the MT directory
    :param panelapp_path: path to the panelapp data dump
    :param config_path: path to the config json
    :param plink_file: pedigree filepath in PLINK format
    :param out_vcf: path to write the VCF out to
    :param tmp_path: temporary checkpoint write path
    """

    # when we checkpoint to a location, we are now reading from it
    # to maintain multiple checkpoints, we rotate between two ext.
    checkpoint_extension = ''
    checkpoint_flip_flop = {'': '2', '2': ''}

    # initiate Hail with defined driver spec.
    init_batch(driver_cores=8, driver_memory='highmem')

    # get the run configuration JSON
    logging.info(f'Reading config dict from "{config_path}"')
    with open(AnyPath(config_path), encoding='utf-8') as handle:
        config_dict = json.load(handle)

    # find the config area specific to hail operations
    hail_config = config_dict.get('filter')

    # read the parsed panelapp data
    logging.info(f'Reading PanelApp data from "{panelapp_path}"')
    with AnyPath(panelapp_path).open() as handle:
        panelapp = json.load(handle)

    # pull green and new genes from the panelapp data
    green_expression, new_expression = green_and_new_from_panelapp(panelapp)

    logging.info(
        f'Starting Hail with reference genome "{hail_config.get("ref_genome")}"'
    )

    # if we already generated the annotated output, load instead
    if not AnyPath(mt_input.rstrip('/') + '/').exists():
        raise Exception(f'Input MatrixTable doesn\'t exist: {mt_input}')

    matrix = hl.read_matrix_table(mt_input)
    logging.debug(
        f'Loaded annotated MT from {mt_input}, size: {matrix.count_rows()}',
    )

    # filter out quality failures
    matrix = filter_on_quality_flags(matrix)

    # running global quality filter steps
    if matrix.count_cols() >= hail_config['min_samples_to_ac_filter']:
        matrix = filter_matrix_by_ac(
            matrix=matrix, ac_threshold=hail_config['ac_threshold']
        )

    matrix = filter_to_well_normalised(matrix)
    matrix = filter_by_ab_ratio(matrix)

    logging.info('checkpointing MT after applying quality filters')
    matrix = matrix.checkpoint(f'{tmp_path}{checkpoint_extension}', overwrite=True)
    checkpoint_extension = checkpoint_flip_flop[checkpoint_extension]
    logging.info(f'Rows remaining: {matrix.count_rows()}')

    # pull annotations into info and update if missing
    logging.info('Pulling VEP annotations into INFO field')
    matrix = extract_annotations(matrix)

    matrix = filter_to_population_rare(matrix=matrix, config=hail_config)
    matrix = split_rows_by_gene_and_filter_to_green(
        matrix=matrix, green_genes=green_expression
    )

    logging.info('checkpointing MT after applying annotation filters')
    matrix = matrix.checkpoint(f'{tmp_path}{checkpoint_extension}', overwrite=True)
    checkpoint_extension = checkpoint_flip_flop[checkpoint_extension]
    matrix = informed_repartition(matrix=matrix)
    logging.info(
        f'Variants remaining after Rare & Green-Gene filter: {matrix.count_rows()}'
    )

    # add Classes to the MT
    # current logic is to apply 1, 2, 3, and 5
    # for cat. 5 (de novo), pre-filter the variants by tx-consequential or C5==1
    logging.info('Applying categories to variant consequences')
    matrix = annotate_category_1(matrix)
    matrix = annotate_category_2(matrix, config=hail_config, new_genes=new_expression)
    matrix = annotate_category_3(matrix, config=hail_config)
    matrix = annotate_category_5(matrix, config=hail_config)
    matrix = annotate_category_4(
        matrix, config=hail_config, plink_family_file=plink_file
    )
    matrix = annotate_category_support(matrix, hail_config)

    # now restrict to only categorised variants
    matrix = filter_to_categorised(matrix)
    matrix = matrix.checkpoint(f'{tmp_path}{checkpoint_extension}', overwrite=True)
    matrix = informed_repartition(matrix=matrix)
    logging.info(f'Variants remaining after Category filter: {matrix.count_rows()}')

    # obtain the massive CSQ string using method stolen from the Broad's Gnomad library
    # also take the single gene_id (from the exploded attribute)
    matrix = matrix.annotate_rows(
        info=matrix.info.annotate(
            CSQ=vep_struct_to_csq(
                matrix.vep, csq_fields=config_dict['variant_object'].get('csq_string')
            ),
            gene_id=matrix.geneIds,
        )
    )

    # write the results to a VCF path
    logging.info(f'Write variants out to "{out_vcf}"')
    write_matrix_to_vcf(
        matrix=matrix,
        output_path=out_vcf,
        additional_header=hail_config.get('csq_header_file'),
    )


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )

    parser = ArgumentParser()
    parser.add_argument(
        '--mt_input',
        required=True,
        help='path to the matrix table to ingest',
    )
    parser.add_argument(
        '--panelapp_path',
        type=str,
        required=True,
        help='bucket path containing panelapp JSON',
    )
    parser.add_argument(
        '--config_path',
        type=str,
        required=False,
        help='If a gene list is being used as a comparison ',
    )
    parser.add_argument(
        '--plink_file', type=str, required=True, help='Family Pedigree in PLINK format'
    )
    parser.add_argument(
        '--out_vcf', type=str, required=True, help='VCF path to export results'
    )
    parser.add_argument(
        '--tmp_path', type=str, required=True, help='Temp path for checkpoints'
    )
    args = parser.parse_args()
    main(
        mt_input=args.mt_input,
        panelapp_path=args.panelapp_path,
        config_path=args.config_path,
        plink_file=args.plink_file,
        out_vcf=args.out_vcf,
        tmp_path=args.tmp_path,
    )
