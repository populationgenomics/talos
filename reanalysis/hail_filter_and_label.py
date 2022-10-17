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
import logging
import sys
from argparse import ArgumentParser
import asyncio
import os

import hail as hl
from peddy import Ped

from cloudpathlib import AnyPath

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import init_batch, output_path, remote_tmpdir, genome_build

# BIG TODO, HACK HERE
def fix_output_path(input: str) -> str:
    return input.replace('/severalgenomes', '/cpg-severalgenomes', 1)
    

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

FIELDS_REQUIRED = {
    'info': [
        ('AC', hl.Int32Expression),
        ('AF', hl.ArrayNumericExpression),
        ('AN', hl.Int32Expression),
    ],
    'splice_ai': [
        ('delta_score', hl.Float32Expression),
        ('splice_consequence', hl.StringExpression),
    ],
    'gnomad_exomes': [
        ('AF', hl.Float64Expression),
        ('AN', hl.Int32Expression),
        ('AC', hl.Int32Expression),
        ('Hom', hl.Int32Expression),
        ('Hemi', hl.Int32Expression),
    ],
    'gnomad_genomes': [
        ('AF', hl.Float64Expression),
        ('AN', hl.Int32Expression),
        ('AC', hl.Int32Expression),
        ('Hom', hl.Int32Expression),
        ('Hemi', hl.Int32Expression),
    ],
    'cadd': [('PHRED', hl.Float32Expression)],
    'dbnsfp': [
        ('REVEL_score', hl.StringExpression),
        ('MutationTaster_pred', hl.StringExpression),
    ],
    'clinvar': [
        ('clinical_significance', hl.StringExpression),
        ('gold_stars', hl.Int32Expression),
    ],
    'geneIds': [],
}

VEP_TX_FIELDS_REQUIRED = [
    ('variant_allele', hl.StringExpression),
    ('consequence_terms', hl.ArrayExpression),
    ('transcript_id', hl.StringExpression),
    ('protein_id', hl.StringExpression),
    ('gene_id', hl.StringExpression),
    ('gene_symbol', hl.StringExpression),
    ('gene_symbol_source', hl.StringExpression),
    ('canonical', hl.Int32Expression),
    ('cdna_start', hl.Int32Expression),
    ('cds_start', hl.Int32Expression),
    ('cds_end', hl.Int32Expression),
    ('biotype', hl.StringExpression),
    ('protein_start', hl.Int32Expression),
    ('protein_end', hl.Int32Expression),
    ('sift_score', hl.Float64Expression),
    ('sift_prediction', hl.StringExpression),
    ('polyphen_score', hl.Float64Expression),
    ('mane_select', hl.StringExpression),
    ('lof', hl.StringExpression),
]


def fields_audit(mt: hl.MatrixTable) -> bool:
    """
    checks that the required fields are all present before continuing
    """
    problems = []
    for field_group, group_types in FIELDS_REQUIRED.items():
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


def vep_audit(mt: hl.MatrixTable) -> bool:
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
        for field, field_type in VEP_TX_FIELDS_REQUIRED:
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


def filter_by_ab_ratio(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    filters HomRef, Het, and HomAlt by appropriate AB ratio bins

    NOTE: This is a broken implementation, as it will replace the
    true genotype calls with missing values. This has implications for
    MOI testing downstream, as the corresponding genotypes in family
    members can be absent, despite being called in the VCF

    This can also cause rows in the resulting VCF to be
    only WT/missing calls, removing all actual variant calls
    :param mt:
    """
    ab = mt.AD[1] / hl.sum(mt.AD)
    return mt.filter_entries(
        (mt.GT.is_hom_ref() & (ab <= 0.15))
        | (mt.GT.is_het() & (ab >= 0.25) & (ab <= 0.75))
        | (mt.GT.is_hom_var() & (ab >= 0.85))
    )


def annotate_category_1(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    applies the boolean Category1 annotation
    semi-rare in Gnomad
    at least one Clinvar star
    contains pathogenic and not conflicting; doesn't contain benign
    :param mt:
    :return: same Matrix, with additional field per variant
    """

    return mt.annotate_rows(
        info=mt.info.annotate(
            categoryboolean1=hl.if_else(
                (mt.info.clinvar_stars > 0)
                & (mt.info.clinvar_sig.lower().contains(PATHOGENIC))
                & ~(mt.info.clinvar_sig.lower().contains(CONFLICTING))
                & ~(mt.info.clinvar_sig.lower().contains(BENIGN)),
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def annotate_category_2(
    mt: hl.MatrixTable, config: dict[str, Any], new_genes: hl.SetExpression
) -> hl.MatrixTable:
    """
    - Gene is new in PanelApp
    - Clinvar contains pathogenic, or
    - Critical protein consequence on at least one transcript
    - High in silico consequence
    :param mt:
    :param config:
    :param new_genes: the new genes in this panelapp content
    :return: same Matrix, with additional field per variant
    """

    critical_consequences = hl.set(config.get('critical_csq'))

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
                    | (
                        (mt.info.clinvar_sig.lower().contains(PATHOGENIC))
                        & ~(mt.info.clinvar_sig.lower().contains(CONFLICTING))
                        & ~(mt.info.clinvar_sig.lower().contains(BENIGN))
                    )
                    | (
                        (mt.info.cadd > config['in_silico']['cadd'])
                        | (mt.info.revel > config['in_silico']['revel'])
                    )
                ),
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def annotate_category_3(mt: hl.MatrixTable, config: dict[str, Any]) -> hl.MatrixTable:
    """
    applies the boolean Category3 flag
    - Critical protein consequence on at least one transcript
    - either predicted NMD or
    - any star Pathogenic or Likely_pathogenic in Clinvar
    :param mt:
    :param config:
    :return:
    """

    critical_consequences = hl.set(config.get('critical_csq'))

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
                    | (
                        (mt.info.clinvar_sig.lower().contains(PATHOGENIC))
                        & ~(mt.info.clinvar_sig.lower().contains(CONFLICTING))
                        & ~(mt.info.clinvar_sig.lower().contains(BENIGN))
                    )
                ),
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def filter_by_consequence(mt: hl.MatrixTable, config: dict[str, Any]) -> hl.MatrixTable:
    """
    - reduce the per-row transcript consequences to a limited group
    - reduce the rows to ones where there are remaining tx consequences
    :param mt:
    :param config: dictionary content relating to hail
    :return: reduced matrix
    """

    # at time of writing this is VEP HIGH + missense_variant
    # update without updating the dictionary content
    high_csq = set(config.get('critical_csq', [])).union(
        set(config.get('additional_consequences', []))
    )

    # overwrite the consequences with an intersection against a limited list
    mt = mt.annotate_rows(
        vep=mt.vep.annotate(
            transcript_consequences=mt.vep.transcript_consequences.filter(
                lambda x: hl.len(hl.set(x.consequence_terms).intersection(high_csq)) > 0
            )
        )
    )

    # filter out rows with no tx consequences left, and no splice cat. assignment
    return mt.filter_rows(
        (hl.len(mt.vep.transcript_consequences) > 0) & (mt.info.categoryboolean5 == 0)
    )


def annotate_category_4(
    mt: hl.MatrixTable, config: dict[str, Any], plink_family_file: str
) -> hl.MatrixTable:
    """
    Category based on de novo MOI, restricted to a group of consequences
    uses the Hail builtin method (very strict)
    :param mt: the whole joint-call MatrixTable
    :param config: all parameters to use when consequence-filtering
    :param plink_family_file: path to a pedigree in PLINK format
    :return: mt with Category4 annotations
    """

    logging.info('Running de novo search')

    de_novo_matrix = filter_by_consequence(mt, config)

    pedigree = hl.Pedigree.read(plink_family_file)

    # avoid consequence filtering twice by calling the de novos in a loop
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


def annotate_category_5(mt: hl.MatrixTable, config: dict[str, Any]) -> hl.MatrixTable:
    """
    :param mt:
    :param config:
    :return: same Matrix, with additional field per variant
    """

    return mt.annotate_rows(
        info=mt.info.annotate(
            categoryboolean5=hl.if_else(
                mt.info.splice_ai_delta >= config['spliceai_full'],
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def annotate_category_support(
    mt: hl.MatrixTable, config: dict[str, Any]
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
    :param mt:
    :param config:
    :return:
    """

    return mt.annotate_rows(
        info=mt.info.annotate(
            categorysupport=hl.if_else(
                (
                    (mt.info.cadd > config['in_silico'].get('cadd'))
                    & (mt.info.revel > config['in_silico'].get('revel'))
                )
                | (
                    (
                        mt.vep.transcript_consequences.any(
                            lambda x: hl.or_else(x.sift_score, MISSING_FLOAT_HI)
                            <= config['in_silico'].get('sift')
                        )
                    )
                    & (
                        mt.vep.transcript_consequences.any(
                            lambda x: hl.or_else(x.polyphen_score, MISSING_FLOAT_LO)
                            >= config['in_silico'].get('polyphen')
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
    mt: hl.MatrixTable, config: dict[str, Any]
) -> hl.MatrixTable:
    """
    run the rare filter, using Gnomad Exomes and Genomes
    :param mt:
    :param config:
    :return:
    """
    # gnomad exomes and genomes below threshold or missing
    # if missing they were previously replaced with 0.0
    return mt.filter_rows(
        (mt.info.gnomad_ex_af < config['af_semi_rare'])
        & (mt.info.gnomad_af < config['af_semi_rare'])
    )


def split_rows_by_gene_and_filter_to_green(
    mt: hl.MatrixTable, green_genes: hl.SetExpression
) -> hl.MatrixTable:
    """
    splits each GeneId onto a new row, then filters any
    rows not annotating a Green PanelApp gene
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
    mt = mt.annotate_rows(
        vep=mt.vep.annotate(
            transcript_consequences=mt.vep.transcript_consequences.filter(
                lambda x: (mt.geneIds == x.gene_id)
                & ((x.biotype == 'protein_coding') | (x.mane_select.contains('NM')))
            )
        )
    )

    return mt


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


def extract_annotations(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    pull out select fields which aren't per-consequence
    store these in INFO (required to be included in VCF export)
    replace with placeholder (least consequential) if empty
    e.g. most tools score 0, but for Sift 1 is least important
    :param mt:
    :return: input matrix with annotations pulled into INFO
    """

    logging.info('Pulling VEP annotations into INFO field')

    return mt.annotate_rows(
        info=mt.info.annotate(
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
            clinvar_sig=hl.or_else(mt.clinvar.clinical_significance, MISSING_STRING),
            clinvar_stars=hl.or_else(mt.clinvar.gold_stars, MISSING_INT),
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
    :param mt:
    :return: input matrix, minus rows without Categories applied
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
    :param mt: the MT to write to file
    """

    # this temp file needs to be in GCP, not local
    # otherwise the batch that generates the file won't be able to read
    additional_cloud_path = fix_output_path(output_path('additional_header.txt', 'tmp'))
    with to_path(additional_cloud_path).open('w') as handle:
        handle.write(
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: '
            'allele|consequence|symbol|gene|feature|mane_select|biotype|exon|hgvsc|'
            'hgvsp|cdna_position|cds_position|protein_position|amino_acids|codons|'
            'allele_num|variant_class|tsl|appris|ccds|ensp|swissprot|trembl|uniparc|'
            'gene_pheno|sift|polyphen|lof|lof_filter|lof_flags">'
        )
    vcf_out = fix_output_path(output_path('hail_categorised.vcf.bgz'))
    logging.info(f'Writing categorised variants out to {vcf_out}')
    hl.export_vcf(mt, vcf_out, append_to_header=additional_cloud_path, tabix=True)


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
    green_genes = set(panel_data.keys()) - {'metadata'}
    logging.info(f'Extracted {len(green_genes)} green genes')
    green_gene_set_expression = hl.literal(green_genes)

    new_genes = {gene for gene in green_genes if panel_data[gene].get('new')}
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


def subselect_mt_to_pedigree(mt: hl.MatrixTable, pedigree: str) -> hl.MatrixTable:
    """
    remove any columns from the MT which are not represented in the Pedigree
    :param mt:
    :param pedigree:
    :return:
    """

    # individual IDs from pedigree
    peddy_ped = Ped(pedigree)
    ped_samples = {individual.sample_id for individual in peddy_ped.samples()}

    # individual IDs from matrix
    matrix_samples = set(mt.s.collect())

    # find overlapping samples
    common_samples = ped_samples.intersection(matrix_samples)

    logging.info(f'Samples in Pedigree: {len(ped_samples)}')
    logging.info(f'Samples in MatrixTable: {len(matrix_samples)}')
    logging.info(f'Common Samples: {len(common_samples)}')

    # full overlap = no filtering
    if common_samples == matrix_samples:
        return mt

    # reduce to those common samples
    mt = mt.filter_cols(hl.literal(common_samples).contains(mt.s))

    logging.info(f'Remaining MatrixTable columns: {mt.count_cols()}')

    return mt


def main(mt_path: str, panelapp: str, config_path: str, plink: str):
    """
    Read the MT from disk
    Do filtering and class annotation
    Export as a VCF
    :param mt_path: path to the MT directory
    :param panelapp: path to the panelapp data dump
    :param config_path: path to the config json
    :param plink: pedigree filepath in PLINK format
    """

    # # initiate Hail with defined driver spec.
    # init_batch(driver_cores=8, driver_memory='highmem')

    asyncio.get_event_loop().run_until_complete(
        hl.init_batch(
            billing_project=get_config()['hail']['billing_project'], 
            remote_tmpdir='hail-az://sevgen002sa/cpg-severalgenomes-hail',
            jar_url="hail-az://hailms02batch/query/jars/1078abac8b8e1c14fe7743aa58bc25118b4108de.jar",
            driver_memory="highmem",
            driver_cores=8
        )
    )


    # checkpoints should be kept independent
    checkpoint_number = 0

    # get the run configuration JSON
    logging.info(f'Reading config dict from "{config_path}"')
    config_dict = read_json_from_path(config_path)

    # get temp suffix from the config (can be None or missing)
    checkpoint_root = fix_output_path(output_path(
        'hail_matrix.mt', config_dict.get('tmp_suffix') or None
    ))

    # find the config area specific to hail operations
    hail_config = config_dict.get('filter')

    # read the parsed panelapp data
    logging.info(f'Reading PanelApp data from "{panelapp}"')
    panelapp = read_json_from_path(panelapp)

    # pull green and new genes from the panelapp data
    green_expression, new_expression = green_and_new_from_panelapp(panelapp)

    logging.info(
        f'Starting Hail with reference genome "{hail_config.get("ref_genome")}"'
    )

    # if we already generated the annotated output, load instead
    if not AnyPath(mt_path.rstrip('/') + '/').exists():
        raise Exception(f'Input MatrixTable doesn\'t exist: {mt_path}')

    mt = hl.read_matrix_table(mt_path)

    if not (fields_audit(mt) and vep_audit(mt)):
        raise Exception('Fields were missing from the input Matrix')

    # subset to currently considered samples
    mt = subselect_mt_to_pedigree(mt, pedigree=plink)

    logging.debug(
        f'Loaded annotated MT from {mt_path}, size: {mt.count_rows()}',
    )

    # filter out quality failures
    mt = filter_on_quality_flags(mt)

    # running global quality filter steps
    mt = filter_matrix_by_ac(mt=mt)
    mt = filter_to_well_normalised(mt)

    # die if there are no variants remaining
    if mt.count_rows() == 0:
        raise Exception('No remaining rows to process!')

    # see method docstring, currently disabled
    # matrix = filter_by_ab_ratio(matrix)

    mt = checkpoint_and_repartition(
        mt,
        checkpoint_root=checkpoint_root,
        checkpoint_num=checkpoint_number,
        extra_logging='after applying quality filters',
    )

    checkpoint_number = checkpoint_number + 1

    mt = extract_annotations(mt)
    mt = filter_to_population_rare(mt=mt, config=hail_config)
    mt = split_rows_by_gene_and_filter_to_green(mt=mt, green_genes=green_expression)

    mt = checkpoint_and_repartition(
        mt,
        checkpoint_root=checkpoint_root,
        checkpoint_num=checkpoint_number,
        extra_logging='after applying Rare & Green-Gene filters',
    )

    checkpoint_number = checkpoint_number + 1

    # add Classes to the MT
    # current logic is to apply 1, 2, 3, and 5, then 4 (de novo)
    # for cat. 4, pre-filter the variants by tx-consequential or C5==1
    logging.info('Applying categories')
    mt = annotate_category_1(mt)
    mt = annotate_category_2(mt, config=hail_config, new_genes=new_expression)
    mt = annotate_category_3(mt, config=hail_config)
    mt = annotate_category_5(mt, config=hail_config)
    mt = annotate_category_4(mt, config=hail_config, plink_family_file=plink)
    mt = annotate_category_support(mt, hail_config)

    mt = filter_to_categorised(mt)
    mt = checkpoint_and_repartition(
        mt,
        checkpoint_root=checkpoint_root,
        checkpoint_num=checkpoint_number,
        extra_logging='after filtering to categorised only',
    )

    # obtain the massive CSQ string using method stolen from the Broad's Gnomad library
    # also take the single gene_id (from the exploded attribute)
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            CSQ=vep_struct_to_csq(
                mt.vep, csq_fields=config_dict['variant_object'].get('csq_string')
            ),
            gene_id=mt.geneIds,
        )
    )

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
    parser.add_argument('--config_path', type=str)
    parser.add_argument('--plink', type=str, required=True, help='Cohort Pedigree')
    args = parser.parse_args()
    main(
        mt_path=args.mt,
        panelapp=args.panelapp,
        config_path=args.config_path,
        plink=args.plink,
    )
