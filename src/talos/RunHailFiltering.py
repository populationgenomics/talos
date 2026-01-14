"""
Read, filter, annotate, classify, and write Genetic data
- read MT
- read PanelApp data through GCP client
- hard-filter (FILTERS, AC)
- extract generic fields
- remove all rows and consequences not relevant to GREEN genes
- extract vep data into CSQ string(s)
- annotate with categories 1, 3, 4, 5, 6, pm5, exomiser, svdb
- remove un-categorised variants
- write as VCF
"""

from argparse import ArgumentParser

from loguru import logger
from mendelbrot.pedigree_parser import PedigreeParser

import hail as hl

from talos.config import config_retrieve
from talos.models import PanelApp
from talos.utils import read_json_from_path

# set some Hail constants
MISSING_INT = hl.int32(0)
MISSING_FLOAT_LO = hl.float64(0.0)
MISSING_STRING = hl.str('missing')
ONE_INT = hl.int32(1)
BENIGN = hl.str('benign')
PATHOGENIC = hl.str('Pathogenic/Likely Pathogenic')
SPLICE_ALTERING = hl.str('splice-altering')
ADDITIONAL_CSQ_DEFAULT = ['missense', 'inframe_deletion', 'inframe_insertion']
CRITICAL_CSQ_DEFAULT = [
    'frameshift',
    'splice_acceptor',
    'splice_donor',
    'start_lost',
    'stop_gained',
    'stop_lost',
    'transcript_ablation',
]
MISSENSE = hl.str('missense')

# decide whether to repartition the data before processing starts
MAX_PARTITIONS = 10000

NUM_PED_COLS = 6


def populate_callset_frequencies(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Annotate the MatrixTable with the callset frequencies
    If these attributes are already present, do nothing
    If AF is the only mising attribute, derive it from AC & AN
    If all are missing, populate using hl.variant_qc
    Returns:
        the original MT, with AC/AF/AN populated
    """

    if all(attribute in mt.info for attribute in ['AC', 'AN', 'AF']):
        logger.info('AC, AN, AF already present, skipping annotation')
        return mt
    if all(attribute in mt.info for attribute in ['AC', 'AN']):
        logger.info('AC, AN present, deriving AF from existing annotations')
        return mt.annotate_rows(
            info=mt.info.annotate(
                AF=hl.map(lambda x: x / mt.info.AN, mt.info.AC),
            ),
        )
    logger.warning('Adding AC/AN/AF annotations to MT based on this callset alone')
    logger.warning('This is unlikely to provide meaningful variant filtering unless this is a huge callset')
    mt = hl.variant_qc(mt, name='variant_qc')
    # per VCF spec, AC and AF should be per alt allele, Hail QC delivers one per allele, including the reference
    # we already mandate that data is normalised/de-duplicated, so we can safely take the last element
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            AC=mt.variant_qc.AC[-1:],
            AN=mt.variant_qc.AN,
            AF=mt.variant_qc.AF[-1:],
        ),
    )
    return mt.drop('variant_qc')


def annotate_clinvarbitration(
    mt: hl.MatrixTable,
    clinvar: str,
    new_genes: hl.SetExpression | None = None,
) -> hl.MatrixTable:
    """
    Don't allow these annotations to be missing
    - Talos has been co-developed with ClinvArbitration, a ClinVar re-summary effort
    - We no longer permit this to be missing (this has slipped in the past, causing odd results)

    See: https://github.com/populationgenomics/ClinvArbitration

    We replace any existing ClinVar annotations with our own version, then identify Pathogenic/Benign variants
    """

    logger.info(f'loading private clinvar annotations from {clinvar}')
    ht = hl.read_table(clinvar)
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            clinvar_significance=hl.or_else(ht[mt.row_key].clinical_significance, MISSING_STRING),
            clinvar_stars=hl.or_else(ht[mt.row_key].gold_stars, MISSING_INT),
            clinvar_allele=hl.or_else(ht[mt.row_key].allele_id, MISSING_INT),
        ),
    )

    # remove all confidently benign
    mt = mt.filter_rows(
        (mt.info.clinvar_significance.lower().contains(BENIGN)) & (mt.info.clinvar_stars > 0),
        keep=False,
    )

    # annotate as either strong or regular, return the result
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            clinvar_talos=hl.if_else(
                mt.info.clinvar_significance == PATHOGENIC,
                ONE_INT,
                MISSING_INT,
            ),
            categorybooleanclinvarplp=hl.if_else(
                (mt.info.clinvar_significance == PATHOGENIC) & (mt.info.clinvar_stars > 0),
                ONE_INT,
                MISSING_INT,
            ),
        ),
    )
    if new_genes is not None:
        mt = mt.annotate_rows(
            info=mt.info.annotate(
                # mark variants that are P/LP, have no stars and in PanelApp "new" genes
                categorybooleanclinvar0starnewgene=hl.if_else(
                    (mt.info.clinvar_significance == PATHOGENIC)
                    & (mt.info.clinvar_stars == 0)
                    & (hl.len(new_genes.intersection(mt.gene_ids)) > 0),
                    ONE_INT,
                    MISSING_INT,
                ),
                # mark variants that are P/LP and have no stars
                categorybooleanclinvar0star=hl.if_else(
                    (mt.info.clinvar_significance == PATHOGENIC) & (mt.info.clinvar_stars == 0),
                    ONE_INT,
                    MISSING_INT,
                ),
            ),
        )
    return mt


def annotate_exomiser(mt: hl.MatrixTable, exomiser: str | None = None, ignored: bool = False) -> hl.MatrixTable:
    """
    Annotate this MT with top hits from Exomiser

    The exomiser table must be indexed on [locus, alleles], with an extra column "proband_details"

    Args:
        mt (): the MatrixTable of all variants
        exomiser (str | None): optional, path to HT of exomiser results
        ignored (bool): if True, don't annotate the MT with the Exomiser data

    Returns:
        The same MatrixTable but with additional annotations
    """
    if not exomiser or ignored:
        logger.info(f'Exomiser not required or requested, skipping annotation (table path: {exomiser})')
        return mt.annotate_rows(info=mt.info.annotate(categorydetailsexomiser=MISSING_STRING))

    logger.info(f'loading exomiser variants from {exomiser}')
    exomiser_ht = hl.read_table(exomiser)
    return mt.annotate_rows(
        info=mt.info.annotate(
            categorydetailsexomiser=hl.or_else(exomiser_ht[mt.row_key].proband_details, MISSING_STRING),
        ),
    )

    logger.info('No exomiser table found, skipping annotation')
    return mt.annotate_rows(info=mt.info.annotate(categorydetailsexomiser=MISSING_STRING))


def annotate_codon_clinvar(mt: hl.MatrixTable, pm5_path: str | None):
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
        pm5_path (str or None): path to a (localised) PM5 annotations Hail Table

    Returns:
        Same MT with an extra category label containing links to all clinvar
        missense variants affecting the same residue as a missense in this
        call-set - shared residue affected on at least one transcript
    """

    if pm5_path is None:
        logger.info(f'PM5 not required or requested, skipping annotation. Table path: {pm5_path}')
        return mt.annotate_rows(info=mt.info.annotate(categorydetailspm5=MISSING_STRING))

    # read in the codon table
    logger.info(f'Reading clinvar alleles by codon from {pm5_path}')
    codon_clinvar = hl.read_table(pm5_path)

    # boom those variants out by consequence
    codon_variants = mt.explode_rows(mt.transcript_consequences).rows()

    # filter for missense, no indels (don't trust VEP)
    codon_variants = codon_variants.filter(
        (hl.len(codon_variants.alleles[0]) == ONE_INT)
        & (hl.len(codon_variants.alleles[1]) == ONE_INT)
        & (codon_variants.transcript_consequences.consequence.contains(MISSENSE)),
    )

    # set the protein residue as an attribute
    codon_variants = codon_variants.annotate(
        residue_affected=hl.str('::').join(
            [
                codon_variants.transcript_consequences.transcript,
                hl.str(codon_variants.transcript_consequences.codon),
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
    return mt.annotate_rows(
        info=mt.info.annotate(
            categorydetailspm5=hl.or_else(codon_variants[mt.row_key].clinvar_variations, MISSING_STRING),
        ),
    )


def annotate_splicevardb(mt: hl.MatrixTable, svdb_path: str | None, ignored: bool = False):
    """
    Takes the locus,ref,alt indexed table of SpliceVarDB variants, and matches
    them up against the variant data
    Annotates with a boolean flag if the variant is Splice-altering according to SVDB

    Args:
        mt (): MT of all variants
        svdb_path (str or None): path to a (localised) PM5 annotations Hail Table
        ignored (bool): if True, don't annotate the MT with the SpliceVarDB data

    Returns:
        Same MT with an extra category label
    """
    if svdb_path is None or ignored:
        logger.info(f'SVDB not required or requested, skipping annotation. (table path: {svdb_path})')
        return mt.annotate_rows(
            info=mt.info.annotate(
                categorybooleansvdb=MISSING_INT,
                svdb_location=MISSING_STRING,
                svdb_method=MISSING_STRING,
                svdb_doi=MISSING_STRING,
            ),
        )

    # read in the codon table
    logger.info(f'Reading SpliceVarDB data from {svdb_path}')
    svdb_ht = hl.read_table(svdb_path)

    # annotate relevant variants with the SVDB results
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            svdb_classification=hl.or_else(svdb_ht[mt.row_key].classification, MISSING_STRING),
            svdb_location=hl.or_else(svdb_ht[mt.row_key].location, MISSING_STRING),
            svdb_method=hl.or_else(svdb_ht[mt.row_key].method, MISSING_STRING),
            svdb_doi=hl.or_else(svdb_ht[mt.row_key].doi, MISSING_STRING),
        ),
    )

    # annotate category if Splice-altering according to SVDB
    return mt.annotate_rows(
        info=mt.info.annotate(
            categorybooleansvdb=hl.if_else(
                mt.info.svdb_classification.lower().contains(SPLICE_ALTERING),
                ONE_INT,
                MISSING_INT,
            ),
        ),
    )


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

    return mt.filter_rows(
        hl.is_missing(mt.filters) | (mt.filters.length() == 0) | (mt.info.clinvar_talos == ONE_INT),
    )


def filter_matrix_by_ac(mt: hl.MatrixTable, af_threshold: float = 0.01, min_ac_to_filter: int = 5) -> hl.MatrixTable:
    """
    Remove variants with AC in joint-call over threshold
    Will never remove variants with 5 or fewer instances
    Also overridden by having a Clinvar Pathogenic anno.

    Expectation here is that the AC and AF values are both arrays, with the first element being the reference

    Args:
        mt (hl.MatrixTable):
        af_threshold (float):
        min_ac_to_filter (int): minimum AC to filter on (no filtering if fewer than X instances in callset)
    Returns:
        MT with all common-in-this-JC variants removed
        (unless overridden by clinvar path)
    """

    return mt.filter_rows(
        ((min_ac_to_filter >= mt.info.AC[0]) | (af_threshold > mt.info.AF[0])) | (mt.info.clinvar_talos == ONE_INT),
    )


def filter_to_population_rare(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    run the rare filter, using Gnomad Exomes and Genomes
    allow clinvar pathogenic to slip through this filter
    """
    # gnomad exomes and genomes below threshold or missing
    # if missing they were previously replaced with 0.0
    # 'semi-rare' as dominant filters will be more strictly filtered later
    rare_af_threshold = config_retrieve(['RunHailFiltering', 'af_semi_rare'])
    return mt.filter_rows(
        (hl.or_else(mt.gnomad.gnomad_AF, MISSING_FLOAT_LO) < rare_af_threshold) | (mt.info.clinvar_talos == ONE_INT),
    )


def remove_variants_outside_gene_roi(mt: hl.MatrixTable, green_genes: hl.SetExpression) -> hl.MatrixTable:
    """
    chunky filtering method - get rid of every variant without at least one green-gene annotation
    does not edit the field itself, that will come later (split_rows_by_gene_and_filter_to_green)

    Args:
        mt ():
        green_genes ():

    Returns:
        the same MT, just reduced
    """

    # filter rows without a green gene (removes empty gene_ids)
    return mt.filter_rows(hl.len(green_genes.intersection(mt.gene_ids)) > 0)


def split_rows_by_gene_and_filter_to_green(mt: hl.MatrixTable, green_genes: hl.SetExpression) -> hl.MatrixTable:
    """
    splits each GeneId onto a new row, then filters any rows not annotating a Green PanelApp gene

    - first explode the matrix, separate gene per row
    - throw away all rows without a green gene
    - on all remaining rows, filter transcript consequences to match _this_ gene

    this is the single most powerful filtering row, effectively leaving just the genes we're interested in

    updated 25/10/2024 - we're retaining snRNA transcripts here, to enable them in

    Args:
        mt ():
        green_genes (): set of all relevant genes
    Returns:
        exploded MatrixTable
    """

    # split each gene onto a separate row, transforms 'gene_ids' field from set to string
    mt = mt.explode_rows(mt.gene_ids)

    # filter rows without a green gene (removes empty gene_ids)
    mt = mt.filter_rows(green_genes.contains(mt.gene_ids))

    # limit the per-row transcript CSQ to those relevant to the single
    # gene now present on each row
    return mt.annotate_rows(
        transcript_consequences=mt.transcript_consequences.filter(
            lambda x: (mt.gene_ids == x.gene_id)
            & ((x.biotype == 'protein_coding') | (x.biotype == 'snRNA') | (x.mane_id.contains('NM'))),
        ),
    )


def annotate_category_alphamissense(mt: hl.MatrixTable) -> hl.MatrixTable:
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
        same variants, categorybooleanalphamissense set to 1 or 0
    """

    return mt.annotate_rows(
        info=mt.info.annotate(
            categorybooleanalphamissense=hl.if_else(
                hl.len(mt.transcript_consequences.filter(lambda x: x.am_class == 'likely_pathogenic')) > 0,
                ONE_INT,
                MISSING_INT,
            ),
        ),
    )


def annotate_category_high_impact(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    applies the boolean Category3 flag
    - Critical protein consequence on at least one transcript

    Args:
        mt (hl.MatrixTable):

    Returns:
        same variants, categorybooleanhighimpact set to 1 or 0
    """

    critical_consequences = hl.set(config_retrieve(['RunHailFiltering', 'critical_csq'], CRITICAL_CSQ_DEFAULT))

    # First check if we have any HIGH consequences
    return mt.annotate_rows(
        info=mt.info.annotate(
            categorybooleanhighimpact=hl.if_else(
                (
                    hl.len(
                        mt.transcript_consequences.filter(
                            lambda x: (
                                hl.len(critical_consequences.intersection(hl.set(x.consequence.split('&')))) > 0
                            ),
                        ),
                    )
                    > 0
                ),
                ONE_INT,
                MISSING_INT,
            ),
        ),
    )


def annotate_category_spliceai(mt: hl.MatrixTable) -> hl.MatrixTable:
    """Label variant with significant spliceAi scores"""
    # skip over MTs where this annotation is absent
    if 'splice_ai' not in mt.row_value:
        return mt.annotate_rows(
            info=mt.info.annotate(
                categorybooleanspliceai=MISSING_INT,
            ),
        )

    # https://github.com/populationgenomics/talos/blob/93bf0455ab233b096a9687aa28d5faf405bac4ef/talos/RunHailFiltering.py#L270C73-L272C105
    return mt.annotate_rows(
        info=mt.info.annotate(
            categorybooleanspliceai=hl.if_else(
                mt.splice_ai.delta_score >= config_retrieve(['RunHailFiltering', 'spliceai']),
                ONE_INT,
                MISSING_INT,
            ),
        ),
    )


def filter_by_consequence(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    - reduce the per-row transcript CSQ to a limited group
    - reduce the rows to ones where there are remaining tx consequences
    """

    # at time of writing this is VEP HIGH + missense_variant
    # update without updating the dictionary content
    critical_consequences = set(config_retrieve(['RunHailFiltering', 'critical_csq'], CRITICAL_CSQ_DEFAULT)) | set(
        config_retrieve(['RunHailFiltering', 'additional_csq'], ADDITIONAL_CSQ_DEFAULT),
    )

    # overwrite the consequences with an intersection against a limited list
    filtered_mt = mt.annotate_rows(
        transcript_consequences=mt.transcript_consequences.filter(
            lambda x: (hl.len(hl.set(x.consequence.split('&')).intersection(critical_consequences)) > 0)
            | (x.biotype == 'snRNA'),
        ),
    )

    # filter out rows with no high impact tx consequences, and no clinvar Pathogenic annotation
    return filtered_mt.filter_rows(
        (hl.len(filtered_mt.transcript_consequences) == 0) & (filtered_mt.info.clinvar_talos == 0),
        keep=False,
    )


def annotate_category_de_novo(
    mt: hl.MatrixTable,
    pedigree_data: PedigreeParser,
    strict_ad: bool = False,
) -> hl.MatrixTable:
    """
    Category based on de novo MOI, restricted to a group of consequences
    The Hail builtin method has limitations around Hemizygous regions

    This was rebuilt by Kyle in the Talkowski lab
    https://github.com/talkowski-lab/variant-interpretation/blob/3a4cdd2e74f0826aa1640473ade03567364602f3/scripts/wes_denovo_full.py

    Both versions heavily rely on probablity likelihoods and AD/Depths to make decisions about whether a call is more
    likely to be noise or a true event.

    In our latest combiner callsets (which may apply at other sites), we don't have PL/AD fields populated for WT calls
    Instead of inserting fake PL and AD values for parents, I'm starting off by not including any of the p value/prior
    testing, using a naive method which checks genotypes and GQ in all family members, AB ratio in child, and permits
    a male proband to have a Hom genotype as a Hemizygous event (looking at you, GATK)

    If this finds the variants we're looking for, great. If we need to reduce noise later, that's fine.
    """

    # modifiable through config
    de_novo_config = config_retrieve(['RunHailFiltering', 'de_novo'])
    min_child_ab: float = de_novo_config.get('min_child_ab', 0.20)
    min_depth: int = de_novo_config.get('min_depth', 5)
    max_depth: int = de_novo_config.get('max_depth', 1000)
    min_proband_gq: int = de_novo_config.get('min_proband_gq', 25)
    min_alt_depth: int = de_novo_config.get('min_alt_depth', 5)

    # this GQ filter will be applied to all samples, not just probands
    # if the dataset is merged from single-sample VCFs, WT samples for a given allele will have no GQ, so this filter
    # will remove all of those samples, removing our ability to detect de novo events in the corresponding families
    min_all_sample_gq: int = de_novo_config.get('min_all_sample_gq', 19)

    # we've seen some pretty weird data formats, so a 'genotype-only' de novo detection option is being added
    # the noise reduction offered by using the more complex filtering strategies have shown great performance internally
    # but if the data you have doesn't allow for those to be used, consider using this flag
    genotype_only: bool = de_novo_config.get('genotype_only', False)

    logger.info('Running de novo search')

    de_novo_matrix = filter_by_consequence(mt)

    # if we use the strict AD filter, we will remove all samples without exactly 2 AD values
    if strict_ad:
        de_novo_matrix = de_novo_matrix.filter_entries(hl.len(de_novo_matrix.AD) == 2)

    # takes the parsed pedigree data, writes it to a temporary file as a Strict 6-column format
    temp_ped_path = 'temp.ped'
    pedigree_data.write_pedigree(output_path=temp_ped_path)
    pedigree = hl.Pedigree.read(temp_ped_path)

    # a series of adjustments for processing sparse data and dragen-igg

    # GT missing -> HomRef if the GQ is above min_all_sample_gq
    de_novo_matrix = de_novo_matrix.annotate_entries(
        GT=hl.if_else(
            # if the GT is assigned, use it
            hl.is_defined(de_novo_matrix.GT),
            de_novo_matrix.GT,
            # if it's a missing value but with a decent GQ, call it a HomRef
            hl.if_else(
                min_all_sample_gq < de_novo_matrix.GQ,
                hl.Call([0, 0]),
                de_novo_matrix.GT,
            ),
        ),
    )

    # pad some AD values if HomRefs are single-element, if AD is missing completely but the var is high quality, insert
    # adjusts for some truly unconventional spec-breaking DRAGEN IGG shenanigans
    de_novo_matrix = de_novo_matrix.annotate_entries(
        AD=hl.if_else(
            # if this is a high quality HomRef
            (de_novo_matrix.GT.is_hom_ref()) & (min_all_sample_gq < de_novo_matrix.GQ),
            hl.if_else(
                # If HomRef is populated
                hl.is_defined(de_novo_matrix.AD),
                hl.if_else(
                    # check for instances of only one AD value, and append a 0 for alt counts
                    (hl.len(de_novo_matrix.AD) == 1) & (de_novo_matrix.AD[0] > 0),
                    de_novo_matrix.AD.append(0),
                    de_novo_matrix.AD,
                ),
                # aiming this at missing "." AD but high quality call, replace with dummy values
                [min_depth + 1, 0],
            ),
            de_novo_matrix.AD,
        ),
    )

    # If DP is 'present' but a missing value, replace it with sum(AD), if DP isn't in the schema at all, insert it
    if 'DP' in de_novo_matrix.entry:
        de_novo_matrix = de_novo_matrix.annotate_entries(
            # if depth is present in the schema and has a real value, use it
            DP=hl.if_else(
                hl.is_defined(de_novo_matrix.DP),
                de_novo_matrix.DP,
                # otherwise use the AD value
                hl.sum(de_novo_matrix.AD),
            ),
        )
    # if depth isn't in the schema at all, insert it
    else:
        de_novo_matrix = de_novo_matrix.annotate_entries(DP=hl.sum(de_novo_matrix.AD))

    # pull out affected members from the pedigree, Hail does not process the phenotype column
    affected_members = hl.literal(pedigree_data.get_affected_member_ids())

    if not genotype_only:
        de_novo_matrix = de_novo_matrix.filter_entries(
            (min_depth > de_novo_matrix.DP)
            | (max_depth < de_novo_matrix.DP)
            | (de_novo_matrix.GT.is_het()) & (de_novo_matrix.AD[1] < (min_child_ab * de_novo_matrix.DP))
            # these tests are aimed exclusively at affected participants
            | ((affected_members.contains(de_novo_matrix.s)) & (min_proband_gq >= de_novo_matrix.GQ)),
            keep=False,
        )

        if min_all_sample_gq:
            logger.info(
                'Applying minimum GQ filter to all samples, this may greatly reduce de Novo event detection if your '
                'dataset was generated from single-sample VCFs. To disable this behaviour, remove '
                '"de_novo.min_all_sample_gq" from the config file.',
            )
            # filter out all samples with GQ below the threshold
            de_novo_matrix = de_novo_matrix.filter_entries(
                (hl.is_defined(de_novo_matrix.GQ)) & (min_all_sample_gq >= de_novo_matrix.GQ),
                keep=False,
            )

    # create a trio matrix (variant rows, trio columns)
    tm = hl.trio_matrix(de_novo_matrix, pedigree, complete_trios=True)

    kid = tm.proband_entry
    dad = tm.father_entry
    mom = tm.mother_entry

    # Updated for hemizygous child calls to not require a call in uninvolved parent
    # call.is_hom_ref is just `all(a == 0 for a in self._alleles)` - no ploidy check
    has_candidate_gt_configuration = (
        (tm.locus.in_autosome_or_par() & kid.GT.is_het() & dad.GT.is_hom_ref() & mom.GT.is_hom_ref())
        | (
            tm.locus.in_x_nonpar()
            & ~tm.is_female
            & ((kid.GT.is_hom_var()) | (kid.GT.is_het()))
            & dad.GT.is_hom_ref()
            & mom.GT.is_hom_ref()
        )
        | (tm.locus.in_x_nonpar() & tm.is_female & kid.GT.is_het() & mom.GT.is_hom_ref() & dad.GT.is_hom_ref())
        | (tm.locus.in_y_nonpar() & ((kid.GT.is_hom_var()) | (kid.GT.is_het())) & dad.GT.is_hom_ref())
        | (tm.locus.in_mito() & kid.GT.is_hom_var() & mom.GT.is_hom_ref())
    )

    # if genotype-only analysis is taking place, the genotypes are all we're using
    if genotype_only:
        tm = tm.annotate_entries(
            de_novo_tested=hl.case().when(has_candidate_gt_configuration, ONE_INT).default(MISSING_INT),
        )

    else:
        # require AD & PL for called variants, just not always for HomRef
        kid_ab = kid.AD[1] / hl.sum(kid.AD)

        # horribly simplified - we don't have the PL or AD for any WTs, so we're really fudging the main parts
        # I've also dropped the requirements for different confidence levels, we're treating Low/Medium/High equally
        tm = tm.annotate_entries(
            de_novo_tested=hl.case()
            .when(~has_candidate_gt_configuration, MISSING_INT)
            .when(min_alt_depth > kid.AD[1], MISSING_INT)
            .when(min_proband_gq > kid.GQ, MISSING_INT)
            .when(kid_ab < min_child_ab, MISSING_INT)
            .default(ONE_INT),
        )

    tm = tm.filter_entries(tm.de_novo_tested == 1)

    # skip most stuff, just retain the keys
    dn_table = tm.entries().select()

    # re-key the table by locus,alleles, removing the sampleID from the compound key
    dn_table = dn_table.key_by(dn_table.locus, dn_table.alleles)

    # we only require the key (locus, alleles) and the sample ID
    # select to remove other fields, then collect per-key into Array of Structs
    dn_table = dn_table.select(dn_table.id).collect_by_key()

    # collect all sample IDs per locus, and squash into a String Array
    # delimit to compress that Array into single Strings
    dn_table = dn_table.annotate(dn_ids=hl.delimit(hl.map(lambda x: x.id, dn_table.values), ','))

    # log the number of variants found this way
    logger.info(f'{dn_table.count()} variants showed de novo inheritance')

    # annotate those values as a flag if relevant, else 'missing'
    return mt.annotate_rows(
        info=mt.info.annotate(categorysampledenovo=hl.or_else(dn_table[mt.row_key].dn_ids, MISSING_STRING)),
    )


def csq_struct_to_string(tx_expr: hl.expr.StructExpression) -> hl.expr.ArrayExpression:
    """
    Taken shamelessly from the gnomad library source code
    Given a CSQ Struct, returns an array of CSQ strings
    (1 per csq in the struct).
    Fields & order correspond to those in `csq_fields`, corresponding to the
    VCF header that is required to interpret the VCF CSQ INFO field.
    Order is flexible & all fields in the default value are supported.
    These fields are formatted in the same way that their VEP CSQ counterparts are.

    Args:
        tx_expr (hl.Struct):
    Returns:
        generates an array of Strings for each CSQ
    """

    def get_csq_from_struct(element: hl.expr.StructExpression) -> hl.expr.StringExpression:
        # Most fields are 1-1, just lowercase
        fields = dict(element)

        # pull the required fields and ordering from config
        csq_fields = config_retrieve(['RunHailFiltering', 'csq_string'])

        return hl.delimit([hl.or_else(hl.str(fields.get(f, '')), '') for f in csq_fields], '|')

    csq = hl.empty_array(hl.tstr)
    csq = csq.extend(
        hl.or_else(tx_expr.map(lambda x: get_csq_from_struct(x)), hl.empty_array(hl.tstr)),
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
        (mt.info.categorybooleanclinvarplp == 1)
        | (mt.info.categorybooleanclinvar0star == 1)
        | (mt.info.categorybooleanclinvar0starnewgene == 1)
        | (mt.info.categorybooleanalphamissense == 1)
        | (mt.info.categorybooleanhighimpact == 1)
        | (mt.info.categorybooleanspliceai == 1)
        | (mt.info.categorysampledenovo != MISSING_STRING)
        | (mt.info.categorydetailspm5 != MISSING_STRING)
        | (mt.info.categorybooleansvdb == 1)
        | (mt.info.categorydetailsexomiser != MISSING_STRING),
    )


def write_matrix_to_vcf(mt: hl.MatrixTable, vcf_out: str):
    """
    write the remaining MatrixTable content to file as a VCF
    generate a custom header containing the CSQ contents which
    were retained during this run

    Args:
        mt (): the whole MatrixTable
        vcf_out (str): where to write the VCF
    """

    # this temp file needs to be in GCP, not local
    # otherwise the batch that generates the file won't be able to read
    header_path = 'additional_header.txt'

    # generate a CSQ string specific to the config file for decoding later
    csq_contents = '|'.join(config_retrieve(['RunHailFiltering', 'csq_string']))

    # write this custom header locally
    with open(header_path, 'w') as handle:
        handle.write(f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: {csq_contents}">')
    logger.info(f'Writing categorised variants out to {vcf_out}')
    hl.export_vcf(mt, vcf_out, append_to_header=header_path, tabix=True)


def green_from_panelapp(panel_data: PanelApp) -> hl.SetExpression:
    """
    Pull all ENSGs from PanelApp data relating to Green Genes

    Args:
        panel_data (PanelApp): the PanelApp object

    Returns:
        a set expression - all Green genes in the analysis
    """

    # take all the green genes, remove the metadata
    green_genes = set(panel_data.genes)
    logger.info(f'Extracted {len(green_genes)} green genes')
    return hl.literal(green_genes)


def new_green_genes_from_panelapp(panel_data: PanelApp) -> hl.SetExpression:
    """
    Pull all ENSGs from PanelApp data labelled as "new" (i.e. PanelDetail.new non-empty)

    Args:
        panel_data (PanelApp): the PanelApp object

    Returns:
        a set expression - all genes marked as new for any panel in the analysis
    """

    new_genes = {ensg for ensg, details in panel_data.genes.items() if getattr(details, 'new', None)}
    logger.info(f'Extracted {len(new_genes)} new genes')
    return hl.literal(new_genes)


def subselect_mt_to_pedigree(mt: hl.MatrixTable, ped_samples: set[str]) -> hl.MatrixTable:
    """Remove any columns from the MT which are not represented in the Pedigree."""

    # individual IDs from matrix
    matrix_samples = set(mt.s.collect())

    # find overlapping samples
    common_samples: set[str] = ped_samples.intersection(matrix_samples)

    logger.info(
        f"""
    Samples in Pedigree: {len(ped_samples)}
    Samples in MatrixTable: {len(matrix_samples)}
    Common Samples: {len(common_samples)}
    """,
    )

    if not common_samples:
        raise ValueError('No samples shared between pedigree and MT')

    # full overlap = no filtering
    if common_samples == matrix_samples:
        return mt

    # reduce to those common samples
    mt = mt.filter_cols(hl.literal(common_samples).contains(mt.s))
    logger.info(f'Remaining MatrixTable columns: {mt.count_cols()}')

    return mt


def generate_a_checkpoint(mt: hl.MatrixTable, checkpoint_path: str) -> hl.MatrixTable:
    """
    wrapper around a few repeated lines which write a checkpoint
    Args:
        mt ():
        checkpoint_path (str): where to write the checkpoint
    """
    logger.info(f'Checkpointing to {checkpoint_path} after filtering out a ton of variants')
    mt = mt.checkpoint(checkpoint_path, overwrite=True)

    # die if there are no variants remaining. Only ever count rows after a checkpoint
    if not (current_rows := mt.count_rows()):
        raise ValueError('No remaining rows to process!')

    logger.info(f'Local checkpoint written, {current_rows} rows remain')
    return mt


def cli_main():
    """
    Read MT, filter, and apply category annotation, export as a VCF
    """

    parser = ArgumentParser()
    parser.add_argument('--input', required=True, help='path to input MT')
    parser.add_argument('--panelapp', help='panelapp JSON')
    parser.add_argument('--pedigree', help='Cohort Pedigree')
    parser.add_argument('--output', help='Where to write the VCF', required=True)
    parser.add_argument('--clinvar', help='HT containing ClinvArbitration annotations', required=True)
    parser.add_argument('--pm5', help='HT containing clinvar PM5 annotations, optional', default=None)
    parser.add_argument('--exomiser', help='HT containing exomiser variant selections, optional', default=None)
    parser.add_argument('--svdb', help='HT containing SpliceVarDB annotations, optional', default=None)
    parser.add_argument('--checkpoint', help='Where/whether to checkpoint, String path', default=None)
    args = parser.parse_args()
    main(
        mt_path=args.input,
        panel_data=args.panelapp,
        pedigree=args.pedigree,
        vcf_out=args.output,
        clinvar=args.clinvar,
        pm5=args.pm5,
        exomiser=args.exomiser,
        svdb=args.svdb,
        checkpoint=args.checkpoint,
    )


def main(  # noqa: PLR0915
    mt_path: str,
    panel_data: str,
    pedigree: str,
    vcf_out: str,
    clinvar: str,
    pm5: str | None = None,
    exomiser: str | None = None,
    svdb: str | None = None,
    checkpoint: str | None = None,
):
    """
    Read MT, filter, and apply category annotation, export as a VCF

    Args:
        mt_path (str): Location of the input MT
        panel_data ():
        pedigree (str): path to a pedigree for this cohort (used in de novo testing)
        vcf_out (str): where to write VCF out
        clinvar (str): path to a ClinVar HT, or unspecified
        pm5 (str): path to a pm5 HT, or unspecified
        exomiser (str): path of an exomiser HT, or unspecified
        svdb (str): path to a SpliceVarDB HT, or unspecified
        checkpoint (str): path to checkpoint data to - serves as checkpoint trigger
    """
    logger.info(
        r"""Welcome To
███████████   █████████   █████          ███████     █████████
█   ███   █  ███     ███   ███         ███     ███  ███     ███
    ███      ███     ███   ███        ███       ███ ███
    ███      ███████████   ███        ███       ███  █████████
    ███      ███     ███   ███        ███       ███         ███
    ███      ███     ███   ███      █  ███     ███  ███     ███
   █████    █████   █████ ███████████    ███████     █████████ """,
    )
    hl.context.init_spark(master='local[*]', default_reference='GRCh38', quiet=True)

    # read the parsed panelapp data
    logger.info(f'Reading PanelApp data from {panel_data!r}')
    panelapp = read_json_from_path(panel_data, return_model=PanelApp)

    # pull green genes from the panelapp data
    green_expression = green_from_panelapp(panelapp)
    # pull genes marked as new in PanelApp
    new_green_genes_expression = new_green_genes_from_panelapp(panelapp)

    # read the pedigree data
    pedigree_data: PedigreeParser = PedigreeParser(pedigree)

    # read the matrix table from a localised directory
    mt = hl.read_matrix_table(mt_path)
    logger.info(f'Loaded annotated MT from {mt_path}, size: {mt.count_rows()}, partitions: {mt.n_partitions()}')

    # Filter out star alleles, not currently capable of handling them
    # Will revisit once our internal experience with DRAGEN-generated variant data improves
    logger.info('Removing any star-allele sites from the dataset, Talos is not currently designed to handle these')
    mt = mt.filter_rows(mt.alleles.contains('*'), keep=False)

    # insert AC/AN/AF if missing
    mt = populate_callset_frequencies(mt)

    # repartition if required - local Hail with finite resources has struggled with some really high (~120k) partitions
    # this creates a local duplicate of the input data with far smaller partition counts, for less processing overhead
    if mt.n_partitions() > MAX_PARTITIONS:
        logger.info('Shrinking partitions way down with an unshuffled repartition')
        mt = mt.repartition(shuffle=False, n_partitions=200)
        if checkpoint:
            logger.info('Trying to write the result locally, might need more space on disk...')
            mt = generate_a_checkpoint(mt, f'{checkpoint}_repartitioned')

    # swap out the default clinvar annotations with private clinvar
    # include a flag for variants that are ClinVar P/LP AND in PanelApp "new" genes
    mt = annotate_clinvarbitration(mt=mt, clinvar=clinvar, new_genes=new_green_genes_expression)

    # remove common-in-gnomad variants (also includes ClinVar annotation)
    mt = filter_to_population_rare(mt=mt)

    # reduce cohort to affected singletons, if the config says so
    if config_retrieve('singletons', False):
        logger.info('Reducing pedigree to affected singletons only')
        pedigree_data.set_participants(pedigree_data.as_singletons())
        pedigree_data.set_participants(pedigree_data.get_affected_members())

    # subset to currently considered samples
    mt = subselect_mt_to_pedigree(mt, ped_samples=pedigree_data.get_all_sample_ids())

    # remove any rows which have no genes of interest
    mt = remove_variants_outside_gene_roi(mt=mt, green_genes=green_expression)

    if checkpoint:
        mt = generate_a_checkpoint(mt, f'{checkpoint}_green_genes')

    # filter out quality failures
    mt = filter_on_quality_flags(mt=mt)

    # filter variants by frequency
    mt = filter_matrix_by_ac(mt=mt)

    # split each gene annotation onto separate rows, filter to green genes (PanelApp ROI)
    mt = split_rows_by_gene_and_filter_to_green(mt=mt, green_genes=green_expression)

    if checkpoint:
        mt = generate_a_checkpoint(mt, f'{checkpoint}_green_and_clean')

    # these categories are ignored early, to prevent an expensive join
    ignored_categories = config_retrieve(['ValidateMOI', 'ignore_categories'], [])

    # annotate this MT with exomiser variants - annotated as MISSING if the table is absent
    mt = annotate_exomiser(mt=mt, exomiser=exomiser, ignored=bool('exomiser' in ignored_categories))

    # if a SVDB data is provided, use that to apply category annotations
    mt = annotate_splicevardb(mt=mt, svdb_path=svdb, ignored=bool('svdb' in ignored_categories))

    # if we ignored both these categories, skip this checkpoint
    if checkpoint and not all(cat in ignored_categories for cat in ['exomiser', 'svdb']):
        mt = generate_a_checkpoint(mt, f'{checkpoint}_green_and_clean_w_external_tables')

    # current logic is to apply 1, 6, 3, then 4 (de novo)
    # 1 was applied earlier during the integration of clinvar data
    # for cat. 4, pre-filter the variants by tx-consequential or C5==1
    logger.info('Applying categories')
    mt = annotate_category_alphamissense(mt=mt)
    mt = annotate_category_high_impact(mt=mt)
    mt = annotate_category_spliceai(mt=mt)

    # insert easy ignore of de novo filtering based on config, to overcome some data format issues
    if any(to_ignore in ignored_categories for to_ignore in ['de_novo', 'denovo', '4']) or config_retrieve(
        'singletons',
        False,
    ):
        logger.info('Skipping de novo annotation, category 4 will not be used during this analysis')
        mt = mt.annotate_rows(info=mt.info.annotate(categorysampledenovo=MISSING_STRING))
    else:
        try:
            # try the standard approach first
            mt = annotate_category_de_novo(mt=mt, pedigree_data=pedigree_data)

        # catch a known error caused by AD fields with a single value
        except hl.utils.java.HailUserError as e:
            logger.error(f'Failed to run de novo annotation, skipping category 4. Error was: {e}')

            # attempt to interpret the error, and re-run with strict AD filtering if it matches
            if 'HailException: array index out of bounds: index=1, length=1' in str(e):
                logger.error(
                    'This error has previously been caused by a variant caller assigning only a single value in the AD '
                    'field for Het calls, where Talos expects two entries; [Ref, Alt].'
                    'If this applies to your dataset, try setting the config option "RunHailFiltering.pad_homref_ad" '
                    'to true, which will add a second value of 0 to all 1-len HomRef AD fields. This will not solve '
                    'situations where the AD field is missing or malformed for Het or Hom calls.'
                    'Talos will re-attempt de novo variant detection, filtering to entries with exactly 2 AD values.',
                )
                mt = annotate_category_de_novo(mt=mt, pedigree_data=pedigree_data, strict_ad=True)
            else:
                mt = mt.annotate_rows(info=mt.info.annotate(categorysampledenovo=MISSING_STRING))

    # if a clinvar-codon table is supplied, use that for PM5
    mt = annotate_codon_clinvar(mt=mt, pm5_path=pm5)

    # remove all variants without positively assigned labels
    mt = filter_to_categorised(mt=mt)

    # obtain the massive CSQ string using method stolen from the Broad's Gnomad library
    # also take the single gene_id (from the exploded attribute)
    # and retrieves the gnomAD annotations
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            **mt.gnomad,
            csq=csq_struct_to_string(mt.transcript_consequences),
            gene_id=mt.gene_ids,
        ),
    )

    write_matrix_to_vcf(mt=mt, vcf_out=vcf_out)


if __name__ == '__main__':
    cli_main()
