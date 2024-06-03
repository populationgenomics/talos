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
from argparse import ArgumentParser
from copy import deepcopy
from datetime import datetime

import backoff
from peds import open_ped

import hail as hl
from hail.utils.java import FatalError

from cpg_utils import to_path
from cpg_utils.config import config_retrieve, get_config, output_path
from cpg_utils.hail_batch import init_batch

from reanalysis.hail_audit import (
    BASE_FIELDS_REQUIRED,
    FIELDS_REQUIRED,
    USELESS_FIELDS,
    VEP_TX_FIELDS_REQUIRED,
    fields_audit,
    vep_audit,
)
from reanalysis.static_values import get_logger
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


def get_clinvar_table(key: str = 'clinvar_decisions') -> str | None:
    """
    try and identify the clinvar table to use
    - try the config specified path
    - fall back to storage:common default path
    - failing that, stick to standard annotations

    Args
        key (str): the key to look for in the config

    Returns:
        a path to a clinvar table, or None
    """

    clinvar_table = get_config()['workflow'].get(key)
    if clinvar_table is not None:
        if to_path(clinvar_table).exists():
            get_logger().info(f'Using clinvar table {clinvar_table}')
            return clinvar_table

    get_logger().info(f'No forced {key} table available, trying default')

    try:
        clinvar_table = to_path(
            os.path.join(
                get_config()['storage']['common']['analysis'],
                'aip_clinvar',
                datetime.now().strftime('%y-%m'),
                f'{key}.ht',
            ),
        )
        # happy path
        if clinvar_table.exists():
            get_logger().info(f'Using clinvar table {clinvar_table}')
            return str(clinvar_table)

        get_logger().info(f'No Clinvar table exists@{clinvar_table}, run the clinvar_runner script')
    except KeyError:
        get_logger().warning('No storage::common::analysis key present')

    return None


def annotate_aip_clinvar(mt: hl.MatrixTable) -> hl.MatrixTable:
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
    Returns:
        The same MatrixTable but with additional annotations
    """

    # if there's private clinvar annotations - use them
    if clinvar := get_clinvar_table():
        # would this replace the standard annotations with missing if there
        # is no private annotation to replace it?
        get_logger().info(f'loading private clinvar annotations from {clinvar}')
        ht = hl.read_table(clinvar)
        mt = mt.annotate_rows(
            info=mt.info.annotate(
                clinvar_significance=hl.or_else(ht[mt.row_key].clinical_significance, MISSING_STRING),
                clinvar_stars=hl.or_else(ht[mt.row_key].gold_stars, MISSING_INT),
                clinvar_allele=hl.or_else(ht[mt.row_key].allele_id, MISSING_INT),
            ),
        )

    # use default annotations
    else:
        get_logger().info('no private annotations, using default contents')

        # do this annotation first, as hail can't string filter against
        # missing contents
        mt = mt.annotate_rows(
            info=mt.info.annotate(
                clinvar_significance=hl.or_else(mt.clinvar.clinical_significance, MISSING_STRING),
                clinvar_stars=hl.or_else(mt.clinvar.gold_stars, MISSING_INT),
                clinvar_allele=hl.or_else(mt.clinvar.allele_id, MISSING_INT),
            ),
        )

    # remove all confidently benign
    mt = mt.filter_rows(
        (mt.info.clinvar_significance.lower().contains(BENIGN)) & (mt.info.clinvar_stars > 0),
        keep=False,
    )

    # annotate as either strong or regular
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            clinvar_aip=hl.if_else(
                (
                    (mt.info.clinvar_significance.lower().contains(PATHOGENIC))
                    & ~(mt.info.clinvar_significance.lower().contains(CONFLICTING))
                ),
                ONE_INT,
                MISSING_INT,
            ),
            clinvar_aip_strong=hl.if_else(
                (
                    (mt.info.clinvar_significance.lower().contains(PATHOGENIC))
                    & ~(mt.info.clinvar_significance.lower().contains(CONFLICTING))
                    & (mt.info.clinvar_stars > 0)
                ),
                ONE_INT,
                MISSING_INT,
            ),
        ),
    )

    return mt


def annotate_codon_clinvar(mt: hl.MatrixTable):
    """
    takes the protein indexed clinvar results and matches up against
    the variant data

    this might be a grossly inefficient method...

    The process is:
    1. load up the hail table of all clinvar annotations indexed by residue
    2. re-shuffle the current variant MT to be similarly indexed on residue
    3. match the datasets across to get residues, clinvar, and corresponding
        loci within this callset
    4. explode that back out to index all the clinvar alleles directly on the
        loci relevant to this callset
    5. use that final table to annotate all relevant clinvar allele IDs as a
        flag, indicating when a variant in this callset creates a residue change
        seen in clinvar

    Note - this might include more munging than necessary, but frankly, I'm not
    smart enough to simplify it. I'll bookmark it for future consideration.

    "Everyone knows that debugging is twice as hard as writing a program in the
    first place. So if you're as clever as you can be when you write it, how will
    you ever debug it?” ― Brian Kernighan

    This matching is universal, i.e. if a variant is the exact position and change
    creating a known pathogenic missense, this method should always find that
    annotation through this lookup. It would be problematic if it didn't. For the
    purposes of PM5, it would be invalid to use the exact same allele as further
    evidence, so exact matches must be filtered out downstream.

    Args:
        mt (): MT of all variants

    Returns:
        Same MT with an extra category label containing links to all clinvar
        missense variants affecting the same residue as a missense in this
        callset - shared residue affected on at least one transcript
    """

    codon_table_path = get_clinvar_table('clinvar_pm5')

    if codon_table_path is None:
        get_logger().info('PM5 table not found, skipping annotation')
        return mt.annotate_rows(info=mt.info.annotate(categorydetailsPM5=MISSING_STRING))

    # read in the codon table
    get_logger().info(f'Reading clinvar alleles by codon from {codon_table_path}')
    codon_clinvar = hl.read_table(str(codon_table_path))

    # boom those variants out by consequence
    codon_variants = mt.explode_rows(mt.vep.transcript_consequences).rows()

    # filter for missense, no indels (don't trust VEP)
    codon_variants = codon_variants.filter(
        (hl.len(codon_variants.alleles[0]) == ONE_INT)
        & (hl.len(codon_variants.alleles[1]) == ONE_INT)
        & (codon_variants.vep.transcript_consequences.consequence_terms.contains('missense_variant')),
    )

    # set the protein residue as an attribute
    codon_variants = codon_variants.annotate(
        residue_affected=hl.str('::').join(
            [
                codon_variants.vep.transcript_consequences.protein_id,
                hl.str(codon_variants.vep.transcript_consequences.protein_start),
            ],
        ),
    )

    # 4. re-key the table on Transcript::Codon
    codon_variants = codon_variants.key_by(codon_variants.residue_affected)

    # 5. extract the position table (protein change linked to all loci)
    codon_variants = codon_variants.select(codon_variants.locus, codon_variants.alleles).collect_by_key(
        name='positions',
    )

    # join the real variant positions with aggregated clinvar
    # 'values' here is the array of all positions
    codon_variants = codon_variants.join(codon_clinvar)

    # explode back out to release the positions
    codon_variants = codon_variants.explode(codon_variants.positions)

    # annotate positions back to normal names (not required?)
    codon_variants = codon_variants.transmute(
        locus=codon_variants.positions.locus,
        alleles=codon_variants.positions.alleles,
    )

    # re-key by locus/allele
    codon_variants = codon_variants.key_by(codon_variants.locus, codon_variants.alleles)

    # aggregate back to position and alleles
    codon_variants = codon_variants.select(codon_variants.clinvar_alleles).collect_by_key(name='clinvar_variations')

    codon_variants = codon_variants.annotate(
        clinvar_variations=hl.str('+').join(
            hl.set(hl.map(lambda x: x.clinvar_alleles, codon_variants.clinvar_variations)),
        ),
    )

    # conditional annotation back into the original MT
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            categorydetailsPM5=hl.or_else(codon_variants[mt.row_key].clinvar_variations, MISSING_STRING),
        ),
    )

    return mt


def filter_on_quality_flags(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    filter MT to rows with 0 quality filters
    note: in Hail, PASS is represented as an empty set

    This is overridden with Clinvar Pathogenic

    Args:
        mt (hl.MatrixTable): all remaining variants
    Returns:
        MT with all filtered variants removed
    """

    return mt.filter_rows(hl.is_missing(mt.filters) | (mt.filters.length() == 0) | (mt.info.clinvar_aip == ONE_INT))


def filter_to_well_normalised(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    single alt per row, no missing Alt

    Args:
        mt (hl.MatrixTable):
    Returns:
        filtered MT
    """

    return mt.filter_rows((hl.len(mt.alleles) == 2) & (mt.alleles[1] != '*'))


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

    get_logger().info('Pulling VEP annotations into INFO field')

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
            splice_ai_csq=hl.or_else(mt.splice_ai.splice_consequence, MISSING_STRING).replace(' ', '_'),
            # we can retain these, but removing completely will lessen the
            # dependence on multiple annotation sources when standing up a
            # new installation
            # cadd=hl.or_else(mt.cadd.PHRED, MISSING_FLOAT_LO),
            # revel=hl.float64(hl.or_else(mt.dbnsfp.REVEL_score, '0.0')),
        ),
    )


def filter_matrix_by_ac(
    mt: hl.MatrixTable,
    ac_threshold: float = 0.01,
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

    return mt.filter_rows(((mt.AC <= 5) | (mt.AC / mt.AN < ac_threshold)) | (mt.info.clinvar_aip == ONE_INT))


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
        ((mt.info.gnomad_ex_af < rare_af_threshold) & (mt.info.gnomad_af < rare_af_threshold))
        | (mt.info.clinvar_aip == ONE_INT),
    )


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
        vep=hl.Struct(transcript_consequences=mt.vep.transcript_consequences, variant_class=mt.vep.variant_class),
    )

    return mt


def split_rows_by_gene_and_filter_to_green(mt: hl.MatrixTable, green_genes: hl.SetExpression) -> hl.MatrixTable:
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
                & ((x.biotype == 'protein_coding') | (x.mane_select.contains('NM'))),
            ),
        ),
    )

    return mt


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
            categoryboolean1=hl.if_else(mt.info.clinvar_aip_strong == ONE_INT, ONE_INT, MISSING_INT),
        ),
    )


def annotate_category_6(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    applies the boolean Category6 flag
    - AlphaMissense likely Pathogenic on at least one transcript
    - Thresholds of am_pathogenicity:
        'Likely benign' if am_pathogenicity < 0.34;
        'Likely pathogenic' if am_pathogenicity > 0.564;
        'ambiguous' otherwise.

    If AM class attribute is missing, skip annotation (default to 0)
    This is run prior to Cat2 annotation - we use C6==True as a replacement
    for the previous approach of [CADD & REVEL predict damaging]

    Args:
        mt (hl.MatrixTable):
    Returns:
        same variants, categoryboolean6 set to 1 or 0
    """

    # focus on the auto-annotated AlphaMissense class
    # allow for the field to be missing
    if 'am_class' not in list(mt.vep.transcript_consequences[0].keys()):
        get_logger().warning('AlphaMissense class not found, skipping annotation')
        return mt.annotate_rows(info=mt.info.annotate(categoryboolean6=MISSING_INT))

    return mt.annotate_rows(
        info=mt.info.annotate(
            categoryboolean6=hl.if_else(
                hl.len(mt.vep.transcript_consequences.filter(lambda x: x.am_class == 'likely_pathogenic')) > 0,
                ONE_INT,
                MISSING_INT,
            ),
        ),
    )


def annotate_category_2(mt: hl.MatrixTable, new_genes: hl.SetExpression | None) -> hl.MatrixTable:
    """
    - Gene is new in PanelApp
    - Clinvar contains pathogenic, or
    - Critical protein consequence on at least one transcript
    - High AlphaMissense score

    Args:
        mt ():
        new_genes (): the new genes in this panelapp content
    Returns:
        same variants, categoryboolean2 set to 1 or 0
    """

    critical_consequences = hl.set(get_config()['filter']['critical_csq'])

    # permit scenario with no new genes
    if new_genes is None:
        return mt.annotate_rows(info=mt.info.annotate(categoryboolean2=MISSING_INT))

    # check for new - if new, allow for in silico, CSQ, or clinvar to confirm
    return mt.annotate_rows(
        info=mt.info.annotate(
            categoryboolean2=hl.if_else(
                (new_genes.contains(mt.geneIds))
                & (
                    (
                        hl.len(
                            mt.vep.transcript_consequences.filter(
                                lambda x: hl.len(critical_consequences.intersection(hl.set(x.consequence_terms))) > 0,
                            ),
                        )
                        > 0
                    )
                    | (mt.info.clinvar_aip == ONE_INT)
                    | (mt.info.categoryboolean6 == ONE_INT)
                ),
                ONE_INT,
                MISSING_INT,
            ),
        ),
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
                            lambda x: (hl.len(critical_consequences.intersection(hl.set(x.consequence_terms))) > 0),
                        ),
                    )
                    > 0
                )
                & (
                    (
                        hl.len(
                            mt.vep.transcript_consequences.filter(
                                lambda x: (hl.len(critical_consequences.intersection(hl.set(x.consequence_terms))) > 0)
                                & ((x.lof == LOFTEE_HC) | (hl.is_missing(x.lof))),
                            ),
                        )
                        > 0
                    )
                    | (mt.info.clinvar_aip == ONE_INT)
                ),
                ONE_INT,
                MISSING_INT,
            ),
        ),
    )


def filter_by_consequence(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    - reduce the per-row transcript CSQ to a limited group
    - reduce the rows to ones where there are remaining tx consequences
    - alternatively except rows with spliceAI delta above threshold

    Args:
        mt ():
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
                lambda x: hl.len(hl.set(x.consequence_terms).intersection(critical_consequences)) > 0,
            ),
        ),
    )

    # filter out rows with no tx consequences left, and no splice cat. assignment
    return filtered_mt.filter_rows(
        (hl.len(filtered_mt.vep.transcript_consequences) == 0) | (filtered_mt.info.categoryboolean5 == 0),
        keep=False,
    )


def annotate_category_4(mt: hl.MatrixTable, ped_file_path: str) -> hl.MatrixTable:
    """
    Category based on de novo MOI, restricted to a group of consequences
    default uses the Hail builtin method (very strict)
    config switch to use the lenient version

    Args:
        mt ():
        ped_file_path (): path to a pedigree in PLINK format

    Returns:
        same variants, categorysample4 either 'missing' or sample IDs
        where de novo inheritance is seen
    """

    get_logger().info('Running de novo search')

    de_novo_matrix = filter_by_consequence(mt)

    pedigree = hl.Pedigree.read(ped_file_path)

    get_logger().info('Updating synthetic PL values for WT calls where missing')

    de_novo_matrix = de_novo_matrix.annotate_entries(
        PL=hl.case()
        .when(~hl.is_missing(de_novo_matrix.PL), de_novo_matrix.PL)
        .when((de_novo_matrix.GT.is_non_ref()) | (hl.is_missing(de_novo_matrix.GQ)), hl.missing('array<int32>'))
        .default([0, de_novo_matrix.GQ, 1000]),
    )

    dn_table = hl.de_novo(
        de_novo_matrix,
        pedigree,
        pop_frequency_prior=de_novo_matrix.info.gnomad_af,
        ignore_in_sample_allele_frequency=True,
        max_parent_ab=get_config()['filter'].get('max_parent_ab', 0.05),
    )

    # re-key the table by locus,alleles, removing the sampleID from the compound key
    dn_table = dn_table.key_by(dn_table.locus, dn_table.alleles)

    # we only require the key (locus, alleles) and the sample ID
    # select to remove other fields, then collect per-key into Array of Structs
    dn_table = dn_table.select(dn_table.id).collect_by_key()

    # collect all sample IDs per locus, and squash into a String Array
    # delimit to compress that Array into single Strings
    dn_table = dn_table.annotate(dn_ids=hl.delimit(hl.map(lambda x: x.id, dn_table.values), ','))

    # log the number of variants found this way
    get_logger().info(f'{dn_table.count()} variants showed de novo inheritance')

    # annotate those values as a flag if relevant, else 'missing'
    return mt.annotate_rows(
        info=mt.info.annotate(**{'categorysample4': hl.or_else(dn_table[mt.row_key].dn_ids, MISSING_STRING)}),
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
            ),
        ),
    )


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

    def get_csq_from_struct(element: hl.expr.StructExpression) -> hl.expr.StringExpression:
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
                'mane_select': element.mane_select,
            },
        )

        # pull the required fields and ordering from config
        csq_fields = get_config()['csq']['csq_string']

        return hl.delimit([hl.or_else(hl.str(fields.get(f, '')), '') for f in csq_fields], '|')

    csq = hl.empty_array(hl.tstr)
    csq = csq.extend(
        hl.or_else(vep_expr['transcript_consequences'].map(lambda x: get_csq_from_struct(x)), hl.empty_array(hl.tstr)),
    )

    # previous consequence filters may make this caution unnecessary
    return hl.or_missing(hl.len(csq) > 0, csq)


def filter_to_categorised(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Filter to rows tagged with a category label

    Args:
        mt ():
    Returns:
        input matrix, minus rows without Categories applied
    """

    return mt.filter_rows(
        (mt.info.categoryboolean1 == 1)
        | (mt.info.categoryboolean6 == 1)
        | (mt.info.categoryboolean2 == 1)
        | (mt.info.categoryboolean3 == 1)
        | (mt.info.categorysample4 != MISSING_STRING)
        | (mt.info.categoryboolean5 == 1)
        | (mt.info.categorydetailsPM5 != MISSING_STRING),
    )


def write_matrix_to_vcf(mt: hl.MatrixTable, vcf_out: str, dataset: str):
    """
    write the remaining MatrixTable content to file as a VCF

    generate a custom header containing the CSQ contents which
    were retained during this run

    Args:
        mt (): the whole MatrixTable
        vcf_out (str): where to write the VCF
        dataset (str): the dataset name
    """

    # this temp file needs to be in GCP, not local
    # otherwise the batch that generates the file won't be able to read
    additional_cloud_path = output_path('additional_header.txt', 'tmp', dataset=dataset)

    # generate a CSQ string specific to the config file for decoding later
    csq_contents = '|'.join(get_config()['csq']['csq_string'])

    # write this custom header locally
    with to_path(additional_cloud_path).open('w') as handle:
        handle.write(f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: {csq_contents}">')
    get_logger().info(f'Writing categorised variants out to {vcf_out}')
    hl.export_vcf(mt, vcf_out, append_to_header=additional_cloud_path, tabix=True)


def green_and_new_from_panelapp(
    panel_genes: dict[str, dict[str, str | list[int]]],
) -> tuple[hl.SetExpression, hl.SetExpression | None]:
    """
    Pull all ENSGs from PanelApp data relating to Green Genes
    Also identify the subset of those genes which relate to NEW in panel

    Args:
        panel_genes (): the 'genes' contents from the panelapp dictionary

    Returns:
        two set expressions - Green, and (Green and New) genes
        can return None if no genes are currently new
    """

    # take all the green genes, remove the metadata
    green_genes = set(panel_genes.keys())
    get_logger().info(f'Extracted {len(green_genes)} green genes')
    green_gene_set_expression = hl.literal(green_genes)

    new_genes = {gene for gene in green_genes if len(panel_genes[gene].get('new', [])) > 0}
    get_logger().info(f'Extracted {len(new_genes)} NEW genes')
    if new_genes:
        return green_gene_set_expression, hl.literal(new_genes)

    return green_gene_set_expression, None


@backoff.on_exception(backoff.expo, exception=FatalError, max_time=60, max_tries=3)
def checkpoint_and_repartition(
    mt: hl.MatrixTable,
    checkpoint_root: str,
    checkpoint_num: int,
    extra_logging: str | None = '',
    prev_mt: str | None = None,
) -> tuple[hl.MatrixTable, str]:
    """
    uses estimated row data size to repartition MT, aiming for a target partition size of ~10MB
    Kat's thread:
    https://discuss.hail.is/t/best-way-to-repartition-heavily-filtered-matrix-tables/2140

    Args:
        mt (): All data
        checkpoint_root (): where to write the checkpoint to
        checkpoint_num (): the checkpoint increment (insert into file path)
        extra_logging (): informative statement to add to logging counts/partitions
    Returns:
        the MT after checkpointing, re-reading, and repartitioning, and new on-disk path
    """
    checkpoint_extended = f'{checkpoint_root}_{checkpoint_num}'
    if (to_path(checkpoint_extended) / '_SUCCESS').exists() and config_retrieve(
        ['workflow', 'reuse_checkpoints'],
        False,
    ):
        get_logger().info(f'Found existing checkpoint at {checkpoint_extended}')
        mt = hl.read_matrix_table(checkpoint_extended)
    else:
        get_logger().info(f'Checkpointing MT to {checkpoint_extended}')
        mt = mt.checkpoint(checkpoint_extended, overwrite=True)

        # once we write and read new data, delete previous
        if prev_mt:
            hl.current_backend().fs.rmtree(prev_mt)

    # estimate partitions; fall back to 1 if low row count
    current_rows = mt.count_rows()
    partitions = current_rows // 200000 or 1

    get_logger().info(f'Re-partitioning {current_rows} into {partitions} partitions {extra_logging}')

    return mt.repartition(n_partitions=partitions, shuffle=True), checkpoint_extended


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
    ped_samples: set[str] = set()
    for family in open_ped(pedigree):
        ped_samples.update({member.id for member in family})

    # individual IDs from matrix
    matrix_samples = set(mt.s.collect())

    # find overlapping samples
    common_samples = ped_samples.intersection(matrix_samples)

    get_logger().info(
        f"""
    Samples in Pedigree: {len(ped_samples)}
    Samples in MatrixTable: {len(matrix_samples)}
    Common Samples: {len(common_samples)}
    """,
    )

    if len(common_samples) == 0:
        raise ValueError('No samples shared between pedigree and MT')

    # full overlap = no filtering
    if common_samples == matrix_samples:
        return mt

    # reduce to those common samples
    mt = mt.filter_cols(hl.literal(common_samples).contains(mt.s))

    get_logger().info(f'Remaining MatrixTable columns: {mt.count_cols()}')

    return mt


def main(
    mt_path: str,
    panelapp_path: str,
    pedigree: str,
    vcf_out: str,
    dataset: str | None = None,
    checkpoint: str | None = None,
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
        checkpoint (str): where to write checkpoints (if at all)
    """

    # never delete this initial data
    initial_mt = deepcopy(mt_path)

    dataset = dataset or get_config()['workflow']['dataset']

    # initiate Hail with defined driver spec.
    init_batch(driver_cores=2, driver_memory='standard')

    # checkpoints should be kept independent
    checkpoint_number = 0

    # get the run configuration JSON
    get_logger().info(f'Reading config dict from {os.getenv("CPG_CONFIG_PATH")}')

    # get temp suffix from the config (can be None or missing)
    # make this checkpoint sequencing-type specific to prevent crossover
    sequencing_type = get_config()['workflow'].get('sequencing_type', 'unknown')
    checkpoint_root = output_path(f'{sequencing_type}_hail_matrix.mt', 'tmp', dataset=dataset)

    # read the parsed panelapp data
    get_logger().info(f'Reading PanelApp data from {panelapp_path!r}')
    panelapp = read_json_from_path(panelapp_path)['genes']  # type: ignore

    # pull green and new genes from the panelapp data
    green_expression, new_expression = green_and_new_from_panelapp(panelapp)

    get_logger().info('Starting Hail with reference genome GRCh38')

    mt = hl.read_matrix_table(mt_path)

    # lookups for required fields all delegated to the hail_audit file
    if not (
        fields_audit(mt=mt, base_fields=BASE_FIELDS_REQUIRED, nested_fields=FIELDS_REQUIRED)  # type: ignore
        and vep_audit(mt=mt, expected_fields=VEP_TX_FIELDS_REQUIRED)
    ):
        mt.describe()
        raise KeyError('Fields were missing from the input Matrix')

    # subset to currently considered samples
    mt = subselect_mt_to_pedigree(mt, pedigree=pedigree)

    get_logger().debug(f'Loaded annotated MT from {mt_path}, size: {mt.count_rows()}')

    # filter out quality failures
    # swap out the default clinvar annotations with private clinvar
    mt = annotate_aip_clinvar(mt=mt)
    mt = filter_on_quality_flags(mt=mt)

    # running global quality filter steps
    mt = filter_to_well_normalised(mt=mt)

    # filter variants by frequency
    mt = extract_annotations(mt=mt)
    mt = filter_matrix_by_ac(mt=mt)
    mt = filter_to_population_rare(mt=mt)

    # shrink the time taken to write checkpoints
    mt = drop_useless_fields(mt=mt)

    if checkpoint:
        # don't delete the starting data
        mt, mt_path = checkpoint_and_repartition(
            mt,
            checkpoint,
            checkpoint_number,
            'after applying quality filters',
            None,
        )
        checkpoint_number += 1

    # die if there are no variants remaining
    assert mt.count_rows(), 'No remaining rows to process!'

    # split genes out to separate rows
    mt = split_rows_by_gene_and_filter_to_green(mt=mt, green_genes=green_expression)

    if checkpoint:
        mt, mt_path = checkpoint_and_repartition(
            mt,
            checkpoint,
            checkpoint_number,
            'after applying Rare & Green-Gene filters',
            prev_mt=mt_path if mt_path != initial_mt else None,
        )
        checkpoint_number += 1

    # add Classes to the MT
    # current logic is to apply 1, 2, 3, and 5, then 4 (de novo)
    # for cat. 4, pre-filter the variants by tx-consequential or C5==1
    get_logger().info('Applying categories')
    mt = annotate_category_1(mt=mt)
    mt = annotate_category_6(mt=mt)
    mt = annotate_category_2(mt=mt, new_genes=new_expression)
    mt = annotate_category_3(mt=mt)
    mt = annotate_category_5(mt=mt)

    # ordering is important - category4 (de novo) makes
    # use of category 5, so it must follow
    mt = annotate_category_4(mt=mt, ped_file_path=pedigree)

    # if a clinvar-codon table is supplied, use that for PM5
    mt = annotate_codon_clinvar(mt=mt)

    mt = filter_to_categorised(mt=mt)

    if checkpoint:
        mt, mt_path = checkpoint_and_repartition(
            mt,
            checkpoint,
            checkpoint_number,
            'after filtering to categorised only',
            prev_mt=mt_path if mt_path != initial_mt else None,
        )
        checkpoint_number += 1

    # obtain the massive CSQ string using method stolen from the Broad's Gnomad library
    # also take the single gene_id (from the exploded attribute)
    mt = mt.annotate_rows(info=mt.info.annotate(CSQ=vep_struct_to_csq(mt.vep), gene_id=mt.geneIds))

    write_matrix_to_vcf(mt=mt, vcf_out=vcf_out, dataset=dataset)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--mt', required=True, help='path to input MT')
    parser.add_argument('--panelapp', type=str, required=True, help='panelapp JSON')
    parser.add_argument('--pedigree', type=str, required=True, help='Cohort Pedigree')
    parser.add_argument('--vcf_out', help='Where to write the VCF', required=True)
    parser.add_argument('--dataset', help='Dataset to write output for')
    parser.add_argument('--checkpoint', help='Path to write checkpoints to', required=False, default=None)
    args = parser.parse_args()
    main(
        mt_path=args.mt,
        panelapp_path=args.panelapp,
        pedigree=args.pedigree,
        vcf_out=args.vcf_out,
        dataset=args.dataset,
        checkpoint=args.checkpoint,
    )
