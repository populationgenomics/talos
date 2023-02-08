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
"""


import os
import logging
import sys
from argparse import ArgumentParser
import asyncio
import os

import hail as hl
from peddy import Ped

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import output_path, init_batch

from reanalysis.hail_audit import (
    fields_audit,
    vep_audit,
    BASE_FIELDS_REQUIRED,
    FIELDS_REQUIRED,
    USELESS_FIELDS,
    VEP_TX_FIELDS_REQUIRED,
)
from reanalysis.utils import read_json_from_path


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


def annotate_aip_clinvar(mt: hl.MatrixTable, clinvar: str) -> hl.MatrixTable:
    """
    instead of making a separate decision about whether the clinvar
    annotation(s) are meaningful during each test, add a single value
    early on. This accomplishes two things:
    1. reduces the logic of each individual test
    2. allows us to use other clinvar consequences for this codon
        (see issue #140)
    3. allows us to replace the clinvar annotation with our private
        annotations (see issue #147)

    Args:
        mt (): the MatrixTable of all variants
        clinvar (str): path to custom table
    Returns:
        The same MatrixTable but with additional annotations
    """

    # if there's private clinvar annotations - use them
    if clinvar != 'absent':
        logging.info(f'loading private clinvar annotations from {clinvar}')
        ht = hl.read_table(clinvar)
        mt = mt.annotate_rows(
            info=mt.info.annotate(
                clinvar_sig=hl.or_else(ht[mt.row_key].rating, MISSING_STRING),
                clinvar_stars=hl.or_else(ht[mt.row_key].stars, MISSING_INT),
                clinvar_allele=hl.or_else(ht[mt.row_key].allele_id, MISSING_INT),
            )
        )

        # remove all confident benign (only confident in this ht)
        mt = mt.filter_rows(mt.info.clinvar_sig.lower().contains(BENIGN), keep=False)

    # use default annotations
    else:
        logging.info(f'no private annotations, using default contents')

        # do this annotation first, as hail can't string filter against
        # missing contents
        mt = mt.annotate_rows(
            info=mt.info.annotate(
                clinvar_sig=hl.or_else(
                    mt.clinvar.clinical_significance, MISSING_STRING
                ),
                clinvar_stars=hl.or_else(mt.clinvar.gold_stars, MISSING_INT),
                clinvar_allele=hl.or_else(mt.clinvar.allele_id, MISSING_INT),
            )
        )

        # remove all confidently benign
        mt = mt.filter_rows(
            (mt.info.clinvar_sig.lower().contains(BENIGN))
            & (mt.info.clinvar_stars > 0),
            keep=False,
        )

    # annotate as either strong or regular
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            clinvar_aip=hl.if_else(
                (
                    (mt.info.clinvar_sig.lower().contains(PATHOGENIC))
                    & ~(mt.info.clinvar_sig.lower().contains(CONFLICTING))
                ),
                ONE_INT,
                MISSING_INT,
            ),
            clinvar_aip_strong=hl.if_else(
                (
                    (mt.info.clinvar_sig.lower().contains(PATHOGENIC))
                    & ~(mt.info.clinvar_sig.lower().contains(CONFLICTING))
                    & (mt.info.clinvar_stars > 0)
                ),
                ONE_INT,
                MISSING_INT,
            ),
        )
    )

    return mt


def filter_matrix_by_ac(
    mt: hl.MatrixTable, ac_threshold: float | None = 0.01
) -> hl.MatrixTable:
    """
    Remove variants with AC in joint-call over threshold
    Will never remove variants with 5 or fewer instances
    Also overridden by having a Clinvar Pathogenic anno.

    Args:
        mt (hl.MatrixTable):
        ac_threshold (float):
    Returns:
        MT with all common-in-this-JC variants removed
        (unless overridden by clinvar path)
    """

    return mt.filter_rows(
        ((mt.AC <= 5) | (mt.AC / mt.AN < ac_threshold))
        | (mt.info.clinvar_aip == ONE_INT)
    )


def filter_on_quality_flags(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    filter MT to rows with 0 quality filters
    note: in Hail, PASS is represented as an empty set

    Args:
        mt (hl.MatrixTable): all remaining variants
    Returns:
        MT with all filtered variants removed
    """

    return mt.filter_rows(mt.filters.length() == 0)


def filter_to_well_normalised(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    single alt per row, no missing Alt

    Args:
        mt (hl.MatrixTable):
    Returns:
        filtered MT
    """

    return mt.filter_rows((hl.len(mt.alleles) == 2) & (mt.alleles[1] != '*'))


def annotate_category_1(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Applies the boolean Category1 annotation
     - clinvar_aip_strong flag, as set in annotate_aip_clinvar
     - represents non-conflicting clinvar pathogenic/likely path

    Args:
        mt ():
    Returns:
        same variants, with categoryboolean1 set to 1 or 0
    """

    return mt.annotate_rows(
        info=mt.info.annotate(
            categoryboolean1=hl.if_else(
                mt.info.clinvar_aip_strong == ONE_INT,
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def annotate_category_2(
    mt: hl.MatrixTable, new_genes: hl.SetExpression
) -> hl.MatrixTable:
    """
    - Gene is new in PanelApp
    - Clinvar contains pathogenic, or
    - Critical protein consequence on at least one transcript
    - High in silico consequence

    Args:
        mt ():
        new_genes (): the new genes in this panelapp content
    Returns:
        same variants, categoryboolean2 set to 1 or 0
    """

    critical_consequences = hl.set(get_config()['filter']['critical_csq'])

    # check for new - if new, allow for in silico, CSQ, or clinvar to confirm
    return mt.annotate_rows(
        info=mt.info.annotate(
            categoryboolean2=hl.if_else(
                (new_genes.contains(mt.geneIds))
                & (
                    (
                        hl.len(
                            mt.vep.transcript_consequences.filter(
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
                    | (mt.info.clinvar_aip == ONE_INT)
                    | (
                        (mt.info.cadd > get_config()['filter']['cadd'])
                        | (mt.info.revel > get_config()['filter']['revel'])
                    )
                ),
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def annotate_category_3(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    applies the boolean Category3 flag
    - Critical protein consequence on at least one transcript
    - either predicted NMD or
    - any star Pathogenic or Likely_pathogenic in Clinvar

    Args:
        mt (hl.MatrixTable):
    Returns:
        same variants, categoryboolean3 set to 1 or 0
    """

    critical_consequences = hl.set(get_config()['filter']['critical_csq'])

    # First check if we have any HIGH consequences
    # then explicitly link the LOFTEE check with HIGH consequences
    # OR allow for a pathogenic ClinVar, any Stars
    return mt.annotate_rows(
        info=mt.info.annotate(
            categoryboolean3=hl.if_else(
                (
                    hl.len(
                        mt.vep.transcript_consequences.filter(
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
                            mt.vep.transcript_consequences.filter(
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
                    | (mt.info.clinvar_aip == ONE_INT)
                ),
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def filter_by_consequence(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    - reduce the per-row transcript CSQ to a limited group
    - reduce the rows to ones where there are remaining tx consequences

    Args:
        mt ():

    Returns:

    """

    # at time of writing this is VEP HIGH + missense_variant
    # update without updating the dictionary content
    critical_consequences = set(get_config()['filter']['critical_csq'])
    additional_consequences = set(get_config()['filter']['additional_csq'])
    critical_consequences.update(additional_consequences)

    # overwrite the consequences with an intersection against a limited list
    filtered_mt = mt.annotate_rows(
        vep=mt.vep.annotate(
            transcript_consequences=mt.vep.transcript_consequences.filter(
                lambda x: hl.len(
                    hl.set(x.consequence_terms).intersection(critical_consequences)
                )
                > 0
            )
        )
    )

    # filter out rows with no tx consequences left, and no splice cat. assignment
    return filtered_mt.filter_rows(
        (hl.len(filtered_mt.vep.transcript_consequences) > 0)
        & (filtered_mt.info.categoryboolean5 == 0)
    )


def annotate_category_4(mt: hl.MatrixTable, plink_family_file: str) -> hl.MatrixTable:
    """
    Category based on de novo MOI, restricted to a group of consequences
    default uses the Hail builtin method (very strict)
    config switch to use the lenient version

    Args:
        mt ():
        plink_family_file (): path to a pedigree in PLINK format

    Returns:
        same variants, categorysample4 either 'missing' or sample IDs
        where de novo inheritance is seen
    """

    logging.info('Running de novo search')

    de_novo_matrix = filter_by_consequence(mt)

    pedigree = hl.Pedigree.read(plink_family_file)

    if get_config()['filter'].get('lenient_de_novo', False):
        logging.info('Inserting synthetic PL values for WT calls')

        # pylint: disable=invalid-unary-operand-type
        de_novo_matrix = de_novo_matrix.annotate_entries(
            PL=hl.case()
            .when(~hl.is_missing(de_novo_matrix.PL), de_novo_matrix.PL)
            .when(
                (de_novo_matrix.GT.is_non_ref()) | (hl.is_missing(de_novo_matrix.GQ)),
                hl.missing('array<int32>'),
            )
            .default([0, de_novo_matrix.GQ, 1000])
        )

    dn_table = hl.de_novo(
        de_novo_matrix,
        pedigree,
        pop_frequency_prior=de_novo_matrix.info.gnomad_af,
        ignore_in_sample_allele_frequency=True,
    )

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
    return mt.annotate_rows(
        info=mt.info.annotate(
            **{
                'categorysample4': hl.or_else(
                    dn_table[mt.row_key].values, MISSING_STRING
                )
            }
        )
    )


def annotate_category_5(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    SpliceAI based category assignment
    Args:
        mt ():

    Returns:
        same variants, categoryboolean5 set to 0 or 1
    """

    return mt.annotate_rows(
        info=mt.info.annotate(
            categoryboolean5=hl.if_else(
                mt.info.splice_ai_delta >= get_config()['filter']['spliceai'],
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def annotate_category_support(mt: hl.MatrixTable) -> hl.MatrixTable:
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

    Args:
        mt ():

    Returns:
        same variants, categorysupport set to 0 or 1
    """

    return mt.annotate_rows(
        info=mt.info.annotate(
            categorysupport=hl.if_else(
                (
                    (mt.info.cadd > get_config()['filter'].get('cadd'))
                    & (mt.info.revel > get_config()['filter'].get('revel'))
                )
                | (
                    (
                        mt.vep.transcript_consequences.any(
                            lambda x: hl.or_else(x.sift_score, MISSING_FLOAT_HI)
                            <= get_config()['filter'].get('sift')
                        )
                    )
                    & (
                        mt.vep.transcript_consequences.any(
                            lambda x: hl.or_else(x.polyphen_score, MISSING_FLOAT_LO)
                            >= get_config()['filter'].get('polyphen')
                        )
                    )
                    & (
                        (mt.info.mutationtaster.contains('D'))
                        | (mt.info.mutationtaster == 'missing')
                    )
                ),
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def filter_to_population_rare(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    run the rare filter, using Gnomad Exomes and Genomes
    allow clinvar pathogenic to slip through this filter
    """
    # gnomad exomes and genomes below threshold or missing
    # if missing they were previously replaced with 0.0
    # 'semi-rare' as dominant filters will be more strictly filtered later
    rare_af_threshold = get_config()['filter']['af_semi_rare']
    return mt.filter_rows(
        (
            (mt.info.gnomad_ex_af < rare_af_threshold)
            & (mt.info.gnomad_af < rare_af_threshold)
        )
        | (mt.info.clinvar_aip == ONE_INT)
    )


def split_rows_by_gene_and_filter_to_green(
    mt: hl.MatrixTable, green_genes: hl.SetExpression
) -> hl.MatrixTable:
    """
    splits each GeneId onto a new row, then filters any
    rows not annotating a Green PanelApp gene

    - first explode the matrix, separate gene per row
    - throw away all rows without a green gene
    - on all remaining rows, filter transript consequences
      to match _this_ gene

    Args:
        mt ():
        green_genes (): set of all relevant genes
    Returns:
        exploded array
    """

    # split each gene onto a separate row
    # transforms 'geneIds' field from set to string
    mt = mt.explode_rows(mt.geneIds)

    # filter rows without a green gene (removes empty geneIds)
    mt = mt.filter_rows(green_genes.contains(mt.geneIds))

    # limit the per-row transcript CSQ to those relevant to the single
    # gene now present on each row
    mt = mt.annotate_rows(
        vep=mt.vep.annotate(
            transcript_consequences=mt.vep.transcript_consequences.filter(
                lambda x: (mt.geneIds == x.gene_id)
                & ((x.biotype == 'protein_coding') | (x.mane_select.contains('NM')))
            )
        )
    )

    return mt


def vep_struct_to_csq(vep_expr: hl.expr.StructExpression) -> hl.expr.ArrayExpression:
    """
    Taken shamelessly from the gnomad library source code
    Given a VEP Struct, returns an array of VEP VCF CSQ strings
    (1 per csq in the struct).
    Fields & order correspond to those in `csq_fields`, corresponding to the
    VCF header that is required to interpret the VCF CSQ INFO field.
    Order is flexible & all fields in the default value are supported.
    These fields are formatted in the same way that their VEP CSQ counterparts are.

    Args:
        vep_expr (hl.Struct):
    Returns:
        generates an array of Strings for each CSQ
    """

    def get_csq_from_struct(
        element: hl.expr.StructExpression,
    ) -> hl.expr.StringExpression:

        # Most fields are 1-1, just lowercase
        fields = dict(element)

        # Add general exceptions
        fields.update(
            {
                'consequence': hl.delimit(element.consequence_terms, delimiter='&'),
                'feature': element.transcript_id,
                'variant_class': vep_expr.variant_class,
                'ensp': element.protein_id,
                'gene': element.gene_id,
                'symbol': element.gene_symbol,
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

        # pull the required fields and ordering from config
        csq_fields = get_config()['csq']['csq_string']

        return hl.delimit(
            [hl.or_else(hl.str(fields.get(f, '')), '') for f in csq_fields], '|'
        )

    csq = hl.empty_array(hl.tstr)
    # pylint: disable=unnecessary-lambda
    csq = csq.extend(
        hl.or_else(
            vep_expr['transcript_consequences'].map(lambda x: get_csq_from_struct(x)),
            hl.empty_array(hl.tstr),
        )
    )

    # previous consequence filters may make this caution unnecessary
    return hl.or_missing(hl.len(csq) > 0, csq)


def extract_annotations(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    pull out select fields which aren't per-consequence
    store these in INFO (required to be included in VCF export)
    replace with placeholder (least consequential) if empty
    e.g. most tools score 0, but for Sift 1 is the least important

    Args:
        mt ():
    Returns:
        Same matrix with re-positioned attributes
    """

    logging.info('Pulling VEP annotations into INFO field')

    return mt.annotate_rows(
        info=mt.info.annotate(
            AC=hl.or_else(mt.AC, MISSING_INT),
            AN=hl.or_else(mt.AN, MISSING_INT),
            gnomad_ex_af=hl.or_else(mt.gnomad_exomes.AF, MISSING_FLOAT_LO),
            gnomad_ex_an=hl.or_else(mt.gnomad_exomes.AN, MISSING_INT),
            gnomad_ex_ac=hl.or_else(mt.gnomad_exomes.AC, MISSING_INT),
            gnomad_ex_hom=hl.or_else(mt.gnomad_exomes.Hom, MISSING_INT),
            gnomad_ex_hemi=hl.or_else(mt.gnomad_exomes.Hemi, MISSING_INT),
            gnomad_af=hl.or_else(mt.gnomad_genomes.AF, MISSING_FLOAT_LO),
            gnomad_an=hl.or_else(mt.gnomad_genomes.AN, MISSING_INT),
            gnomad_ac=hl.or_else(mt.gnomad_genomes.AC, MISSING_INT),
            gnomad_hom=hl.or_else(mt.gnomad_genomes.Hom, MISSING_INT),
            gnomad_hemi=hl.or_else(mt.gnomad_genomes.Hemi, MISSING_INT),
            splice_ai_delta=hl.or_else(mt.splice_ai.delta_score, MISSING_FLOAT_LO),
            splice_ai_csq=hl.or_else(
                mt.splice_ai.splice_consequence, MISSING_STRING
            ).replace(' ', '_'),
            cadd=hl.or_else(mt.cadd.PHRED, MISSING_FLOAT_LO),
            # these next 3 are per-transcript, with ';' to delimit
            # pulling these annotations into INFO with ';' to separate
            # will break INFO parsing for most tools
            revel=hl.float64(hl.or_else(mt.dbnsfp.REVEL_score, '0.0')),
            mutationtaster=hl.or_else(
                mt.dbnsfp.MutationTaster_pred, MISSING_STRING
            ).replace(';', ','),
        )
    )


def filter_to_categorised(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Filter to rows tagged with a class

    Args:
        mt ():
    Returns:
        input matrix, minus rows without Categories applied
    """

    return mt.filter_rows(
        (mt.info.categoryboolean1 == 1)
        | (mt.info.categoryboolean2 == 1)
        | (mt.info.categoryboolean3 == 1)
        | (mt.info.categorysample4 != 'missing')
        | (mt.info.categoryboolean5 == 1)
        | (mt.info.categorysupport == 1)
    )


def write_matrix_to_vcf(mt: hl.MatrixTable):
    """
    write the remaining MatrixTable content to file as a VCF

    Args:
        mt (): the whole MatrixTable
    Returns:
        path to write MT out to
    """

    # this temp file needs to be in GCP, not local
    # otherwise the batch that generates the file won't be able to read
    additional_cloud_path = output_path('additional_header.txt', 'tmp')
    with to_path(additional_cloud_path).open('w') as handle:
        handle.write(
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: '
            'allele|consequence|symbol|gene|feature|mane_select|biotype|exon|hgvsc|'
            'hgvsp|cdna_position|cds_position|protein_position|amino_acids|codons|'
            'allele_num|variant_class|tsl|appris|ccds|ensp|swissprot|trembl|uniparc|'
            'gene_pheno|sift|polyphen|lof|lof_filter|lof_flags">'
        )
    vcf_out = output_path('hail_categorised.vcf.bgz', 'analysis')
    logging.info(f'Writing categorised variants out to {vcf_out}')
    hl.export_vcf(mt, vcf_out, append_to_header=additional_cloud_path, tabix=True)


def green_and_new_from_panelapp(
    panel_genes: dict[str, dict[str, str]]
) -> tuple[hl.SetExpression, hl.SetExpression]:
    """
    Pull all ENSGs from PanelApp data relating to Green Genes
    Also identify the subset of those genes which relate to NEW in panel

    Args:
        panel_genes (): the 'genes' contents from the panelapp dictionary

    Returns:
        two set expressions - Green, and (Green and New) genes
    """

    # take all the green genes, remove the metadata
    green_genes = set(panel_genes.keys())
    logging.info(f'Extracted {len(green_genes)} green genes')
    green_gene_set_expression = hl.literal(green_genes)

    new_genes = {
        gene for gene in green_genes if len(panel_genes[gene].get('new', [])) > 0
    }
    logging.info(f'Extracted {len(new_genes)} NEW genes')
    new_gene_set_expression = hl.literal(new_genes)

    return green_gene_set_expression, new_gene_set_expression


def checkpoint_and_repartition(
    mt: hl.MatrixTable,
    checkpoint_root: str,
    checkpoint_num: int,
    extra_logging: str | None = '',
) -> hl.MatrixTable:
    """
    uses estimated row data size to repartition MT
    aiming for a target partition size of ~10MB
    Kat's thread:
    https://discuss.hail.is/t/best-way-to-repartition-heavily-filtered-matrix-tables/2140

    Args:
        mt (): All data
        checkpoint_root (): where to write the checkpoint to
        checkpoint_num (): the checkpoint increment (insert into file path)
        extra_logging (): informative statement to add to logging counts/partitions
    Returns:
        the MT after checkpointing, re-reading, and repartitioning
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


def subselect_mt_to_pedigree(mt: hl.MatrixTable, pedigree: str) -> hl.MatrixTable:
    """
    remove any columns from the MT which are not represented in the Pedigree
    _not_ done: remove variants without calls amongst the current samples
    that's probably more processing than benefit - will be skipped later

    Args:
        mt ():
        pedigree ():
    Returns:
        MatrixTable, possibly with fewer columns (samples) present
    """

    # individual IDs from pedigree
    peddy_ped = Ped(pedigree)
    ped_samples = {individual.sample_id for individual in peddy_ped.samples()}

    # individual IDs from matrix
    matrix_samples = set(mt.s.collect())

    # find overlapping samples
    common_samples = ped_samples.intersection(matrix_samples)

    logging.info(
        f"""
    Samples in Pedigree: {len(ped_samples)}
    Samples in MatrixTable: {len(matrix_samples)}
    Common Samples: {len(common_samples)}
    """
    )

    if len(common_samples) == 0:
        raise ValueError('No samples shared between pedigree and MT')

    # full overlap = no filtering
    if common_samples == matrix_samples:
        return mt

    # reduce to those common samples
    mt = mt.filter_cols(hl.literal(common_samples).contains(mt.s))

    logging.info(f'Remaining MatrixTable columns: {mt.count_cols()}')

    return mt


def drop_useless_fields(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Remove fields from the MT which are irrelevant to the analysis
    Write times for the checkpoints are a substantial portion of the runtime
    Aim to reduce that by removing the amount of written data

    These fields would not be exported in the VCF anyway, so no downstream
    impacts caused by removal prior to that write

    Args:
        mt ():
    Returns:

    """

    # drop the useless top-level fields
    mt = mt.drop(*[field for field in USELESS_FIELDS if field in mt.row_value])

    # now drop most VEP fields
    mt = mt.annotate_rows(
        vep=hl.Struct(
            transcript_consequences=mt.vep.transcript_consequences,
            variant_class=mt.vep.variant_class,
        )
    )

    return mt


def main(mt_path: str, panelapp: str, plink: str, clinvar: str):
    """
    Read MT, filter, and apply category annotation
    Export as a VCF
    """

    # # initiate Hail with defined driver spec.
    init_batch(driver_cores=8, driver_memory='highmem')

    # checkpoints should be kept independent
    checkpoint_number = 0

    # get the run configuration JSON
    logging.info(f'Reading config dict from {os.getenv("CPG_CONFIG_PATH")}')

    # get temp suffix from the config (can be None or missing)
    checkpoint_root = output_path('hail_matrix.mt', 'tmp')

    # read the parsed panelapp data
    logging.info(f'Reading PanelApp data from {panelapp!r}')
    panelapp = read_json_from_path(panelapp)['genes']

    # pull green and new genes from the panelapp data
    green_expression, new_expression = green_and_new_from_panelapp(panelapp)

    logging.info('Starting Hail with reference genome GRCh38')

    # if we already generated the annotated output, load instead
    if not to_path(mt_path.rstrip('/') + '/').exists():
        raise FileExistsError(f'Input MatrixTable doesn\'t exist: {mt_path}')

    mt = hl.read_matrix_table(mt_path)

    # lookups for required fields all delegated to the hail_audit file
    if not (
        fields_audit(
            mt=mt, base_fields=BASE_FIELDS_REQUIRED, nested_fields=FIELDS_REQUIRED
        )
        and vep_audit(mt=mt, expected_fields=VEP_TX_FIELDS_REQUIRED)
    ):
        raise KeyError('Fields were missing from the input Matrix')

    # subset to currently considered samples
    mt = subselect_mt_to_pedigree(mt, pedigree=plink)

    logging.debug(
        f'Loaded annotated MT from {mt_path}, size: {mt.count_rows()}',
    )

    # filter out quality failures
    mt = filter_on_quality_flags(mt=mt)

    # running global quality filter steps
    mt = filter_to_well_normalised(mt=mt)

    # shrink the time taken to write checkpoints
    mt = drop_useless_fields(mt=mt)

    # mt = checkpoint_and_repartition(
    #     mt=mt,
    #     checkpoint_root=checkpoint_root,
    #     checkpoint_num=checkpoint_number,
    #     extra_logging='after applying quality filters',
    # )

    # die if there are no variants remaining
    if mt.count_rows() == 0:
        raise ValueError('No remaining rows to process!')

    # checkpoint_number = checkpoint_number + 1

    # swap out the default clinvar annotations with private clinvar
    mt = annotate_aip_clinvar(mt=mt, clinvar=clinvar)
    mt = extract_annotations(mt=mt)

    # filter variants by frequency
    mt = filter_matrix_by_ac(mt=mt)
    mt = filter_to_population_rare(mt=mt)

    # split genes out to separate rows
    mt = split_rows_by_gene_and_filter_to_green(mt=mt, green_genes=green_expression)

    # mt = checkpoint_and_repartition(
    #     mt=mt,
    #     checkpoint_root=checkpoint_root,
    #     checkpoint_num=checkpoint_number,
    #     extra_logging='after applying Rare & Green-Gene filters',
    # )

    # checkpoint_number = checkpoint_number + 1

    # add Classes to the MT
    # current logic is to apply 1, 2, 3, and 5, then 4 (de novo)
    # for cat. 4, pre-filter the variants by tx-consequential or C5==1
    logging.info('Applying categories')
    mt = annotate_category_1(mt=mt)
    mt = annotate_category_2(mt=mt, new_genes=new_expression)
    mt = annotate_category_3(mt=mt)
    mt = annotate_category_5(mt=mt)

    # cat. 4 can run in 2 modes - config contains a switch
    mt = annotate_category_4(mt=mt, plink_family_file=plink)
    mt = annotate_category_support(mt=mt)

    mt = filter_to_categorised(mt=mt)
    # mt = checkpoint_and_repartition(
    #     mt=mt,
    #     checkpoint_root=checkpoint_root,
    #     checkpoint_num=checkpoint_number,
    #     extra_logging='after filtering to categorised only',
    # )

    # obtain the massive CSQ string using method stolen from the Broad's Gnomad library
    # also take the single gene_id (from the exploded attribute)
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            CSQ=vep_struct_to_csq(mt.vep),
            gene_id=mt.geneIds,
        )
    )

    logging.info('Writing VCF')
    write_matrix_to_vcf(mt=mt)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )

    parser = ArgumentParser()
    parser.add_argument('--mt', required=True, help='path to input MT')
    parser.add_argument('--panelapp', type=str, required=True, help='panelapp JSON')
    parser.add_argument('--plink', type=str, required=True, help='Cohort Pedigree')
    parser.add_argument(
        '--clinvar', type=str, default='absent', help='Custom Clinvar Summary HT'
    )
    args = parser.parse_args()
    main(
        mt_path=args.mt, panelapp=args.panelapp, plink=args.plink, clinvar=args.clinvar
    )
