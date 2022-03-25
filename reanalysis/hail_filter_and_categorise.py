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

from typing import Any, Dict, Optional, Tuple
from itertools import permutations
import json
import logging
import sys

import click
import hail as hl

from cloudpathlib import AnyPath


# set some Hail constants
MISSING_STRING = hl.str('missing')
MISSING_INT = hl.int32(0)
ONE_INT = hl.int32(1)
MISSING_FLOAT_LO = hl.float64(0.0)
MISSING_FLOAT_HI = hl.float64(1.0)

CONFLICTING = hl.str('conflicting')
LOFTEE_HC = hl.str('HC')
PATHOGENIC = hl.str('pathogenic')


def annotate_category_1(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    applies the Category1 annotation (1 or 0) as appropriate
    semi-rare in Gnomad
    at least one Clinvar star
    either Pathogenic or Likely_pathogenic in Clinvar

    Didn't handle 'Pathogenic/Likely_pathogenic'
    Changing to 'contains pathogenic and not conflicting'
    :param matrix:
    :return: same Matrix, with additional field per variant
    """

    return matrix.annotate_rows(
        info=matrix.info.annotate(
            Category1=hl.if_else(
                (matrix.info.clinvar_stars > 0)
                & (matrix.info.clinvar_sig.lower().contains(PATHOGENIC))
                & ~(matrix.info.clinvar_sig.lower().contains(CONFLICTING)),
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def annotate_category_2(
    matrix: hl.MatrixTable, config: Dict[str, Any], new_genes: hl.SetExpression
) -> hl.MatrixTable:
    """
    - Gene is new in PanelApp
    - Rare in Gnomad, and
    - Clinvar, or
    - Critical protein consequence on at least one transcript
    - High in silico consequence

    New update! this is now restricted to the NEW genes only
    This means that these are now confident Category2, and we
    only have a MOI test remaining

    :param matrix:
    :param config:
    :param new_genes: the new genes in this panelapp content
    :return: same Matrix, with additional field per variant
    """

    critical_consequences = hl.set(config.get('critical_csq'))

    return matrix.annotate_rows(
        info=matrix.info.annotate(
            Category2=hl.if_else(
                (new_genes.contains(matrix.geneIds))
                & (
                    (
                        matrix.vep.transcript_consequences.any(
                            lambda x: hl.len(
                                critical_consequences.intersection(
                                    hl.set(x.consequence_terms)
                                )
                            )
                            > 0
                        )
                    )
                    | (matrix.info.clinvar_sig.lower().contains(PATHOGENIC))
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
    matrix: hl.MatrixTable, config: Dict[str, Any]
) -> hl.MatrixTable:
    """
    applies the Category3 flag where appropriate
    - Critical protein consequence on at least one transcript
    - rare in Gnomad
    - either predicted NMD or
    - any star Pathogenic or Likely_pathogenic in Clinvar
    :param matrix:
    :param config:
    :return:
    """

    critical_consequences = hl.set(config.get('critical_csq'))

    return matrix.annotate_rows(
        info=matrix.info.annotate(
            Category3=hl.if_else(
                (
                    matrix.vep.transcript_consequences.any(
                        lambda x: hl.len(
                            critical_consequences.intersection(
                                hl.set(x.consequence_terms)
                            )
                        )
                        > 0
                    )
                )
                & (
                    (
                        matrix.vep.transcript_consequences.any(
                            lambda x: (x.lof == LOFTEE_HC) | (hl.is_missing(x.lof))
                        )
                    )
                    | (matrix.info.clinvar_sig.lower().contains(PATHOGENIC))
                ),
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def annotate_category_4(
    matrix: hl.MatrixTable, config: Dict[str, Any]
) -> hl.MatrixTable:
    """
    Class based on in silico annotations
    - rare in Gnomad, and
    - CADD & REVEL above threshold (switched to consensus), or
    - Massive cross-tool consensus
    - polyphen and sift are evaluated per-consequence
    :param matrix:
    :param config:
    :return:
    """

    return matrix.annotate_rows(
        info=matrix.info.annotate(
            Category4=hl.if_else(
                (
                    (matrix.info.cadd > config.get('cadd'))
                    & (matrix.info.revel > config.get('revel'))
                )
                | (
                    (
                        matrix.vep.transcript_consequences.any(
                            lambda x: hl.or_else(x.sift_score, MISSING_FLOAT_HI)
                            <= config.get('sift')
                        )
                    )
                    & (
                        matrix.vep.transcript_consequences.any(
                            lambda x: hl.or_else(x.polyphen_score, MISSING_FLOAT_LO)
                            >= config.get('polyphen')
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


def annotate_category_4_only(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    applies a flag to all variants with only Category4 flags applied
    this becomes relevant when we are looking at variants eligible for
    compound het analysis

    :param matrix:
    """

    return matrix.annotate_rows(
        category_4_only=hl.if_else(
            (matrix.info.Category1 == 0)
            & (matrix.info.Category2 == 0)
            & (matrix.info.Category3 == 0)
            & (matrix.info.Category4 == 1),
            ONE_INT,
            MISSING_INT,
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
        class_4_only=0
    )

    transform into simplified 1-10-GC-G
    drop the class_4_only attribute
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


def extract_comp_het_details(matrix: hl.MatrixTable):
    """
    takes the matrix table, and finds compound-hets per sample
    based on the gene name only

    :param matrix:
    """

    # set a new group of values as the key, so that we can collect on them easily
    ch_matrix = matrix.key_rows_by(matrix.locus, matrix.alleles, matrix.class_4_only)
    ch_matrix = ch_matrix.annotate_cols(
        hets=hl.agg.group_by(
            ch_matrix.info.gene_id,
            hl.agg.filter(ch_matrix.GT.is_het(), hl.agg.collect(ch_matrix.row_key)),
        )
    )

    # extract those possible compound het pairs out as a non-Hail structure
    compound_hets = {}

    # iterate over the hail table rows
    # find all variant pair permutations which aren't both class 4
    for row in ch_matrix.select_cols('hets').col.collect():

        # prepare a summary dict for this sample
        sample_dict = {}

        # iterate over all the `gene: [var1, var2]` structures
        for gene, variants in dict(row.hets).items():

            # assess each possible variant pairing
            for var1, var2 in permutations(variants, 2):

                # skip if both are class 4 only - not valuable pairing
                if var1.class_4_only == 1 and var2.class_4_only == 1:
                    continue

                # pair the string transformation
                sample_dict.setdefault(gene, []).append(
                    [transform_variant_string(var1), transform_variant_string(var2)]
                )

        # if we found comp hets, add the content for this sample
        if len(sample_dict) > 0:
            compound_hets[row.s] = sample_dict

    return compound_hets


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
    failure conditions applied are dependent on whether VQSR was run
    if VQSR - allow for filters to be empty, PASS, or VQSR with a PASS category assigned
    if not - require the variant filters to be PASS
    :param matrix_data:
    :param vqsr_run: if True, we
    :return:
    """
    pass_literal = hl.literal('PASS')
    pass_set = hl.literal({pass_literal})
    vqsr_literal = hl.literal('VQSR')
    vqsr_set = hl.literal({vqsr_literal})

    if vqsr_run:

        # hard filter for quality; assuming data is well normalised in pipeline
        matrix_data = matrix_data.filter_rows(
            (
                (matrix_data.filters.length() == 0)
                | (matrix_data.filters == pass_set)
                | (
                    (matrix_data.filters == vqsr_set)
                    & (matrix_data.info.AS_FilterStatus)
                )
            )
        )

    # otherwise strictly enforce FILTERS==PASS
    else:
        matrix_data = matrix_data.filter_rows(matrix_data.filters == pass_set)

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


def annotate_using_vep(matrix_data: hl.MatrixTable) -> hl.MatrixTable:
    """
    runs VEP annotation on the MT
    :param matrix_data:
    :return:
    """

    # throw in a repartition here (annotate even chunks in parallel)
    matrix_data = matrix_data.repartition(150, shuffle=False)

    # VEP 105 annotations
    return hl.vep(matrix_data, config='file:///vep_data/vep-gcloud.json')


def filter_rows_for_rare(
    matrix: hl.MatrixTable, config: Dict[str, Any]
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


def filter_benign(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    filter out benign variants, where clinvar is confident
    not worth running this filter separately
    :param matrix:
    """
    # remove all clinvar benign, decent level of support
    benign = hl.str('benign')
    return matrix.filter_rows(
        (matrix.info.clinvar_sig.lower().contains(benign))
        & (matrix.info.clinvar_stars > 0),
        keep=False,
    )


def filter_to_green_genes_and_split(
    matrix: hl.MatrixTable, green_genes: hl.SetExpression
) -> hl.MatrixTable:
    """
    reduces geneIds set to green only, then splits
    :param matrix:
    :param green_genes:
    """

    # remove any rows with no genic annotation at all
    # do this here as set intersections with Missing will fail
    matrix = matrix.filter_rows(hl.is_missing(matrix.geneIds), keep=False)

    # replace the default list of green IDs with a reduced set
    matrix = matrix.annotate_rows(geneIds=green_genes.intersection(matrix.geneIds))

    # split to form a separate row for each green gene
    # this transforms the 'geneIds' field from a set to a string
    return matrix.explode_rows(matrix.geneIds)


def filter_by_consequence(
    matrix: hl.MatrixTable, config: Dict[str, Any]
) -> hl.MatrixTable:
    """
    - reduce the per-row transcript consequences to those specific to the geneIds
    - reduce the rows to ones where there are remaining tx consequences

    :param matrix:
    :param config: dictionary content relating to hail
    :return: reduced matrix
    """

    # identify consequences to discard from the config
    useless_csq = hl.set(config['useless_csq'])

    # reduce consequences to overlap with per-variant green geneIDs (pre-filtered)
    # added another condition to state that the tx biotype needs to be protein_coding,
    # unless the row also has an attached MANE transcript
    # consider an extra allowance for strong Appris transcripts
    matrix = matrix.annotate_rows(
        vep=matrix.vep.annotate(
            transcript_consequences=matrix.vep.transcript_consequences.filter(
                lambda x: (matrix.geneIds == x.gene_id)
                & (hl.len(hl.set(x.consequence_terms).difference(useless_csq)) > 0)
                & ((x.biotype == 'protein_coding') | (x.mane_select.contains('NM')))
            )
        )
    )

    # filter out all rows with no remaining consequences
    return matrix.filter_rows(hl.len(matrix.vep.transcript_consequences) > 0)


def vep_struct_to_csq(
    vep_expr: hl.expr.StructExpression, csq_fields: str
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
    :param csq_fields: | delimited fields to include in the CSQ (in that order)
    :return: The corresponding CSQ strings
    """
    _csq_fields = [f.lower() for f in csq_fields.split('|')]

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
            [hl.or_else(hl.str(fields.get(f, '')), '') for f in _csq_fields], '|'
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
        | (matrix.info.Category4 == 1)
    )


def write_matrix_to_vcf(matrix: hl.MatrixTable, output_path: str):
    """
    write the remaining MatrixTable content to file as a VCF
    :param matrix:
    :param output_path: where to write
    """
    hl.export_vcf(
        matrix,
        output_path,
        tabix=True,
    )


def green_and_new_from_panelapp(
    panel_data: Dict[str, Dict[str, str]]
) -> Tuple[hl.SetExpression, hl.SetExpression]:
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


@click.command()
@click.option('--mt', 'mt_path', help='path to the matrix table to ingest')
@click.option('--pap', 'panelapp_path', help='bucket path containing panelapp JSON')
@click.option('--config', 'config_path', help='path to a config dict')
@click.option('--output', 'out_vcf', help='VCF path to export results')
@click.option(
    '--mt_out',
    help='path to export annotated MT to',
    required=False,
    default=None,
)
def main(
    mt_path: str,
    panelapp_path: str,
    config_path: str,
    out_vcf: str,
    mt_out: Optional[str] = None,
):
    """
    Read the MT from disk
    Do filtering and class annotation
    Export as a VCF

    :param mt_path: path to the MT directory
    :param panelapp_path: path to the panelapp data dump
    :param config_path: path to the config json
    :param out_vcf: path to write the VCF out to
    :param mt_out:
    """

    # get the run configuration JSON
    logging.info(f'Reading config dict from "{config_path}"')
    with open(AnyPath(config_path), encoding='utf-8') as handle:
        config_dict = json.load(handle)

    # find the config area specific to hail operations
    hail_config = config_dict.get('filter')

    # read the parsed panelapp data
    logging.info(f'Reading PanelApp data from "{panelapp_path}"')
    with open(AnyPath(panelapp_path), encoding='utf-8') as handle:
        panelapp = json.load(handle)

    # pull green and new genes from the panelapp data
    green_expression, new_expression = green_and_new_from_panelapp(panelapp)

    logging.info(
        f'Starting Hail with reference genome "{hail_config.get("ref_genome")}"'
    )
    # initiate Hail with the specified reference
    hl.init(default_reference=hail_config.get('ref_genome'), quiet=True)

    # if we already generated the annotated output, load instead
    if AnyPath(mt_out.rstrip('/') + '/').exists():
        logging.info(
            f'Loading annotated MT from "{mt_out}"',
        )
        matrix = hl.read_matrix_table(mt_out)
        logging.debug(
            f'Loaded annotated MT, size: {matrix.count_rows()}',
        )

    else:
        logging.info(f'Loading MT from "{mt_path}"')
        matrix = hl.read_matrix_table(mt_path)
        logging.debug(f'Loaded pre-annotation MT, size: {matrix.count_rows()}')

        # hard filter entries in the MT prior to annotation
        logging.info('Hard filtering variants')
        matrix = filter_matrix_by_ac(matrix_data=matrix, config=hail_config)
        matrix = filter_matrix_by_variant_attributes(matrix_data=matrix)
        logging.info(f'Post-Hard filtering size: {matrix.count_rows()}')

        # re-annotate using VEP
        logging.info('Annotating variants')
        matrix = annotate_using_vep(matrix_data=matrix)

        # if a path is provided, dump the MT
        if mt_out is not None:
            matrix.write(mt_out, overwrite=True)

    # pull annotations into info and update if missing
    logging.info('Pulling VEP annotations into INFO field')
    matrix = extract_annotations(matrix)

    # filter on row annotations
    logging.info('Filtering Variant rows')
    matrix = filter_rows_for_rare(matrix=matrix, config=hail_config)
    logging.info(f'Variants remaining after Rare filter: {matrix.count_rows()}')
    matrix = filter_benign(matrix=matrix)
    logging.info(f'Variants remaining after Benign filter: {matrix.count_rows()}')
    matrix = filter_to_green_genes_and_split(
        matrix=matrix, green_genes=green_expression
    )
    logging.info(f'Variants remaining after Green-Gene filter: {matrix.count_rows()}')
    matrix = filter_by_consequence(matrix=matrix, config=hail_config)
    logging.info(f'Variants remaining after Consequence filter: {matrix.count_rows()}')

    # choose some logical way of repartitioning
    logging.info('Repartition to 50 fragments following Gene ID filter')
    matrix = matrix.repartition(50, shuffle=True)

    # add Classes to the MT
    logging.info('Applying classes to variant consequences')
    matrix = annotate_category_1(matrix)
    matrix = annotate_category_2(matrix, hail_config, new_expression)
    matrix = annotate_category_3(matrix, hail_config)
    matrix = annotate_category_4(matrix, hail_config['in_silico'])

    # filter to class-annotated only prior to export
    logging.info('Filter variants to leave only classified')
    matrix = filter_to_categorised(matrix)
    logging.info(f'Variants remaining after Category filter: {matrix.count_rows()}')

    # add an additional annotation, if the variant is Category4 only
    matrix = annotate_category_4_only(matrix)

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

    # parse out the compound het details (after pulling gene_id above)
    comp_het_details = extract_comp_het_details(matrix=matrix)

    # transform the vcf output path into a json path
    out_json = f'{out_vcf.split(".", maxsplit=1)[0]}.json'

    # and write the comp-het JSON file
    serialised_obj = json.dumps(comp_het_details, indent=True, default=str)
    AnyPath(out_json).write_text(serialised_obj)

    logging.info('comp-het data written to cloud')

    # write the results to a VCF path
    logging.info(f'Write variants out to "{out_vcf}"')
    write_matrix_to_vcf(matrix=matrix, output_path=out_vcf)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%M-%d %H:%M:%S',
        stream=sys.stderr,
    )
    main()  # pylint: disable=E1120
