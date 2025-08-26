"""
classes and methods shared across reanalysis components

HTTPX requests are backoff-wrapped using tenacity
https://tenacity.readthedocs.io/en/latest/
"""

import json
import pathlib
import re
import sqlite3
import statistics
import string
import zoneinfo
from collections import defaultdict
from datetime import datetime
from itertools import chain, combinations, combinations_with_replacement, islice
from random import choices
from typing import TYPE_CHECKING, Any

import httpx
from cloudpathlib.anypath import to_anypath
from loguru import logger
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential_jitter

from talos.config import config_retrieve
from talos.models import (
    VARIANT_MODELS,
    Coordinates,
    ResultData,
    SmallVariant,
    StructuralVariant,
    lift_up_model_version,
    translate_category,
)
from talos.pedigree_parser import PedigreeParser
from talos.static_values import get_granular_date

if TYPE_CHECKING:
    import cyvcf2


HOMREF: int = 0
HETALT: int = 1
UNKNOWN: int = 2
HOMALT: int = 3

TWO = 2
THREE = 3

# in cyVCF2, these ints represent HOMREF, and UNKNOWN
BAD_GENOTYPES: set[int] = {HOMREF, UNKNOWN}
PHASE_SET_DEFAULT = -2147483648
NON_HOM_CHROM = ['X', 'Y', 'MT', 'M']
X_CHROMOSOME = {'X'}

TIMEZONE = zoneinfo.ZoneInfo('Australia/Brisbane')
TODAY = datetime.now(tz=TIMEZONE).strftime('%Y-%m-%d_%H:%M')

# most lenient to most conservative
# usage = if we have two MOIs for the same gene, take the broadest
ORDERED_MOIS = ['Mono_And_Biallelic', 'Monoallelic', 'Hemi_Mono_In_Female', 'Hemi_Bi_In_Female', 'Biallelic']
IRRELEVANT_MOI = {'unknown', 'other'}

DATE_RE = re.compile(r'\d{4}-\d{2}-\d{2}')

# this just saves some typing
MEMBER_LOOKUP_DICT = {'0': None}

DOI_URL = 'https://doi.org/'

# we've noted some instances where we failed the whole process due to failure to parse phase data
# don't fail, just suggest that someone raises an issue on GitHub, but only print this message once
PHASE_BROKEN: bool = False


# template SQL for interacting with the history/state DB
INSERT_VAR_QUERY = 'INSERT INTO variants (contig, position, reference, alternate) VALUES (?, ?, ?, ?);'
INSERT_CAT_QUERY = 'INSERT INTO categories (var_id, sample_id, category, date) VALUES (?, ?, ?, ?);'
INSERT_PARTNER_QUERY = 'INSERT INTO partners (var_id, sample_id, partners) VALUES (?, ?, ?);'
INSERT_VARSTARS_QUERY = 'INSERT INTO var_stars (var_id, sample_id, clinvar_stars) VALUES (?, ?, ?);'

# update the clinvar stars if we see a new value
UPDATE_CLINVAR_QUERY = 'UPDATE var_stars SET clinvar_stars = ? WHERE var_id = ? AND sample_id = ?;'

# inner join on variants and categories, left join on partners (optional). Basically, get everything for a sample
QUERY_SAMPLE_ALL = """
    SELECT
        variants.var_id,
        variants.contig || '-' || variants.position || '-' || variants.reference || '-' || variants.alternate as var_key,
        categories.category,
        categories.date,
        var_stars.clinvar_stars,
        partners.partners
    FROM
        (variants INNER JOIN categories ON variants.var_id = categories.var_id)
        LEFT JOIN partners ON variants.var_id = partners.var_id AND categories.sample_id = partners.sample_id
        LEFT JOIN var_stars ON variants.var_id = var_stars.var_id AND categories.sample_id = var_stars.sample_id
    WHERE
        categories.sample_id = "{}";
"""
QUERY_ALL_VAR_IDS = (
    "SELECT var_id, contig || '-' || position || '-' || reference || '-' || alternate as var_key FROM variants;"
)
QUERY_ALL_VAR_STARS = """
SELECT 
    variants.var_id, 
    variants.contig || '-' || variants.position || '-' || variants.reference || '-' || variants.alternate as var_key,
    var_stars.sample_id, 
    var_stars.first_pheno_match 
FROM variants INNER JOIN var_stars ON variants.var_id = var_stars.var_id;
"""
INSERT_VARSTARS_PHENO = 'INSERT INTO var_stars (var_id, sample_id, first_pheno_match) VALUES (?, ?, ?);'
UPDATE_VARSTARS_PHENO = 'UPDATE var_stars SET first_pheno_match = ? WHERE var_id = ? AND sample_id = ?;'


def parse_mane_json_to_dict(mane_json: str) -> dict:
    """
    Read the MANE JSON and filter it to the relevant fields

    Returns:
        a dictionary of {Symbol: ID}
    """

    json_dict = read_json_from_path(mane_json)

    return {entry['symbol']: entry['ensg'] for entry in json_dict.values()}


def get_random_string(length: int = 6) -> str:
    """
    get a random string of a pre-determined leng`th

    Returns:
        A random string comprised of upper-case letters and numbers
    """
    return ''.join(choices(string.ascii_uppercase + string.digits, k=length))


def chunks(iterable, chunk_size):
    """
    Yield successive n-sized chunks from an iterable

    Args:
        iterable (): any iterable - tuple, str, list, set
        chunk_size (): size of intervals to return

    Returns:
        intervals of requested size across the collection
    """

    if isinstance(iterable, set):
        iterable = list(iterable)

    for i in range(0, len(iterable), chunk_size):
        yield iterable[i : (i + chunk_size)]


def generator_chunks(generator, size):
    """
    Iterates across a generator, returning specifically sized chunks

    Args:
        generator (): any generator or method implementing yield
        size (): size of iterator to return

    Returns:
        a subset of the generator results
    """
    iterator = iter(generator)
    for first in iterator:
        yield list(chain([first], islice(iterator, size - 1)))


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential_jitter(initial=1, max=5, exp_base=2),
    retry=retry_if_exception_type(
        (
            httpx.ReadTimeout,
            httpx.ConnectError,
            httpx.TooManyRedirects,
            httpx.RequestError,
        ),
    ),
    reraise=True,
)
def get_json_response(url):
    """
    takes a request URL, checks for healthy response, returns the JSON
    For this purpose we only expect a dictionary return
    List use-case (activities endpoint) no longer supported

    for backoff handling I'm using the tenacity library

    Args:
        url (str): URL to retrieve JSON format data from

    Returns:
        the JSON response from the endpoint
    """
    response = httpx.get(url, headers={'Accept': 'application/json'}, timeout=60, follow_redirects=True)
    if response.is_success:
        return response.json()
    raise ValueError('The JSON response could not be parsed successfully')


def get_phase_data(samples: list[str], var: 'cyvcf2.Variant') -> dict[str, dict[int, str]]:
    """
    read phase data from this variant

    we tolerate a variety of failures, the worst case scenario is that we have no phase data
    PS - PhaseSet
    PID/PGT - PhaseID/PhaseGenotype

    No other phase types are currently supported

    Args:
        samples (list[str]): all samples in the VCF
        var (cyvcf2.Variant):
    """

    global PHASE_BROKEN

    phased_dict: dict[str, dict[int, str]] = defaultdict(dict)

    # check we have relevant attributes in the variant format fields
    if not any(x in var.FORMAT for x in ['PS', 'PGT', 'PID']):
        return dict(phased_dict)

    # first set the numpy.ndarray to be a list of ints
    # then zip against ordered sample IDs
    # this might need to store the exact genotype too
    # i.e. 0|1 and 1|0 can be in the same phase-set
    # but are un-phased variants
    try:
        if 'PS' in var.FORMAT:
            for sample, phase, genotype in zip(samples, map(int, var.format('PS')), var.genotypes, strict=True):
                # cyvcf2.Variant holds two ints, and a bool for biallelic calls
                # but only one int and a bool for hemi
                if len(genotype) == THREE:
                    allele_1, allele_2, phased = genotype
                    gt = f'{allele_1}|{allele_2}'
                elif len(genotype) == TWO:
                    allele_1, phased = genotype
                    gt = f'{allele_1}'
                else:
                    raise ValueError(f'Unexpected genotype length: {len(genotype)}: {genotype}')

                if not phased:
                    continue

                # phase set is a number
                if phase != PHASE_SET_DEFAULT:
                    phased_dict[sample][phase] = gt

        elif all(attr_name in var.FORMAT for attr_name in ['PGT', 'PID']):
            logger.info('Failed to find PS phase attributes')
            try:
                # retry using PGT & PID
                for sample, phase_gt, phase_id in zip(
                    samples,
                    var.format('PGT'),
                    var.format('PID'),
                    strict=True,
                ):
                    if phase_gt != '.' and phase_id != '.':
                        phased_dict[sample][phase_id] = phase_gt

            except KeyError as ke2:
                logger.info('Failed to determine phase information using PID and PGT')
                raise ke2
        elif not PHASE_BROKEN:
            logger.info('Found no PS phase attributes (known formats are PS, PGT/PID)')
            PHASE_BROKEN = True

    except (KeyError, ValueError):
        if not PHASE_BROKEN:
            logger.info('Failed to correctly parse known phase attributes using existing methods')
            logger.info(
                'Please post an issue on the Talos GitHub Repo with the VCF FORMAT lines and descriptions',
            )
        PHASE_BROKEN = True

    return dict(phased_dict)


def organise_exomiser(
    info_dict: dict[str, Any],
    rank_threshold: int | None = None,
):
    """
    method dedicated to handling the new exomiser annotations

    If present, these are in a condensed format of PROBAND_RANK_MOI, delimited by "::" e.g.
    "PROBAND1_29_AR::PROBAND2_26_AR::PROBAND3_23_AR"

    When parsing this, we optionally filter to only the top-n ranked results

    Args:
        info_dict ():
        rank_threshold (int | None): if present, only retain results above this rank
    """

    info_dict['categorysampleexomiser'] = []

    # this becomes a dict of dicts - Family, MOI, rank
    info_dict['exomiser'] = defaultdict(dict)

    # if completely absent, the 'samples' category annotation is an empty set
    if 'categorydetailsexomiser' not in info_dict:
        return

    # pop off the exomiser details
    exomiser_details = info_dict.pop('categorydetailsexomiser')

    if exomiser_details == 'missing':
        return

    # split the string into a list of strings, iterate over the list
    for each_exomiser in exomiser_details.split('::'):
        # split the string into its component parts
        sam_id, rank, moi = each_exomiser.split('_')

        # swap that rank to an int
        rank_int: int = int(rank)

        # optionally, filter by a threshold rank
        if rank_threshold and rank_int > rank_threshold:
            continue

        info_dict['exomiser'][sam_id][moi] = rank_int

    # keep the sample-type category for use in "is this variant labelled for this sample" checks
    # it's not exact - we need the pedigree for this to work at all
    # and there is the potential for catching family members without the variant
    info_dict['categorysampleexomiser'] = list(info_dict['exomiser'].keys())


def polish_exomiser_results(results: ResultData) -> None:
    """
    now that we have all the per-participant events, we pare back any exomiser content in the variant info
    field to be exclusively those relevant to each sample

    Args:
        results (ResultData): the results object to be updated
    """
    # iterate over the samples in the results part of the data model
    for sample, content in results.results.items():
        # for each variant
        for var in content.variants:
            # if this sample had exomiser ratings
            if 'exomiser' in var.categories and isinstance(var.var_data.info['exomiser'], defaultdict):
                # pull out the exomiser section specific to this sample
                exomiser_section = var.var_data.info['exomiser'][sample]
                # shove them in a list of strings
                var.var_data.info['exomiser'] = [f'{moi}:{rank}' for moi, rank in exomiser_section.items()]
            else:
                # if it wasn't categorised for this sample, empty list
                var.var_data.info['exomiser'] = []


def organise_pm5(info_dict: dict[str, Any]):
    """
    method dedicated to handling the new pm5 annotations

    e.g. categorydetailsPM5=27037::Pathogenic::1+27048::Pathogenic::1;
    1. break into component allele data

    Returns:
        None, updates self. attributes
    """

    if 'categorydetailspm5' not in info_dict:
        return

    pm5_content = info_dict.pop('categorydetailspm5')

    # nothing to do here
    if pm5_content == 'missing':
        info_dict['categorybooleanpm5'] = 0
        return

    # current clinvar annotation, if any
    current_clinvar = str(info_dict.get('clinvar_allele', 'not_this'))

    # instantiate a dict to store csq-matched results
    pm5_data = {}

    # break the strings into a set
    pm5_strings = set(pm5_content.split('+'))
    for clinvar_entry in pm5_strings:
        # fragment each entry
        allele_id, stars = clinvar_entry.split('::')

        # never consider the exact match, pm5 is always separate
        if allele_id == current_clinvar:
            continue

        # if non-self, add to the dict
        pm5_data[allele_id] = stars

    # case where no non-self alleles were found
    # assigning False and not-assigning are equivalent, just return
    if pm5_data:
        # set boolean category and specific data
        info_dict['categorybooleanpm5'] = 1
        info_dict['pm5_data'] = pm5_data

    else:
        info_dict['categorybooleanpm5'] = 0


def organise_svdb_doi(info_dict: dict[str, Any]):
    """
    method dedicated to handling the SV DB DOI records
    edits in place

    Args:
        info_dict ():
    """
    if 'svdb_doi' not in info_dict:
        return

    # pop off the value
    doi_value = info_dict.pop('svdb_doi')

    if doi_value == 'missing':
        info_dict['svdb_doi'] = []
        return

    # split the value
    doi_urls = []
    for doi in doi_value.split(','):
        doi_urls.append(DOI_URL + doi)
    info_dict['svdb_doi'] = doi_urls


def create_small_variant(
    var: 'cyvcf2.Variant',
    samples: list[str],
):
    """
    takes a small variant and creates a Model from it

    Args:
        var ():
        samples ():
    """

    coordinates = Coordinates(chrom=var.CHROM.replace('chr', ''), pos=var.POS, ref=var.REF, alt=var.ALT[0])
    info: dict[str, Any] = {x.lower(): y for x, y in var.INFO} | {'var_link': coordinates.string_format}

    het_samples, hom_samples = get_non_ref_samples(variant=var, samples=samples)

    # organise PM5
    organise_pm5(info)

    # organise SVDB DOIs
    organise_svdb_doi(info)

    # organise the exomiser data, if present. By default, only retain the top 2 ranked results
    organise_exomiser(
        info,
        rank_threshold=config_retrieve(
            ['ValidateMOI', 'exomiser_rank_threshold'],
            2,
        ),
    )

    # optionally - ignore some categories from this analysis
    ignored_categories = config_retrieve(['ValidateMOI', 'ignore_categories'], [])

    # set the class attributes - skipping over categories we've chosen to ignore
    boolean_categories = [
        key
        for key in info
        if key.startswith('categoryboolean') and key.replace('categoryboolean', '') not in ignored_categories
    ]
    sample_categories = [
        key
        for key in info
        if key.startswith('categorysample') and key.replace('categorysample', '') not in ignored_categories
    ]

    # the categories to be treated as support-only for this runtime - make it a set
    support_categories = set(config_retrieve(['ValidateMOI', 'support_categories'], []))

    # overwrite with true booleans
    for cat in boolean_categories:
        info[cat] = info.get(cat, 0) == 1

    # sample categories are a list of strings or 'missing'
    for sam_cat in sample_categories:
        if isinstance(info[sam_cat], str):
            info[sam_cat] = info[sam_cat].split(',') if info[sam_cat] != 'missing' else []
        elif isinstance(info[sam_cat], list):
            info[sam_cat] = info[sam_cat]
        elif isinstance(info[sam_cat], set):
            info[sam_cat] = list(info[sam_cat])

    phased = get_phase_data(samples, var)

    # only keep these where the sample has a variant - the majority of samples have empty data, and we don't use it
    # if we require depths/ratios/etc. for WT samples, revisit this
    # hopefully a solution to the memory explosion in large cohorts
    variant_samples = het_samples | hom_samples
    depths: dict[str, int] = {
        k: v for k, v in zip(samples, map(int, var.gt_depths), strict=True) if k in variant_samples
    }
    alt_depths: dict[str, int] = {
        k: v for k, v in zip(samples, map(int, var.gt_alt_depths), strict=True) if k in variant_samples
    }
    ab_ratios: dict[str, float] = {
        k: v for k, v in zip(samples, map(float, var.gt_alt_freqs), strict=True) if k in variant_samples
    }
    transcript_consequences = extract_csq(csq_contents=info.pop('csq', ''))

    return SmallVariant(
        coordinates=coordinates,
        info=info,
        het_samples=het_samples,
        hom_samples=hom_samples,
        boolean_categories=boolean_categories,
        sample_categories=sample_categories,
        ignored_categories=ignored_categories,
        support_categories=support_categories,
        phased=phased,
        alt_depths=alt_depths,
        depths=depths,
        ab_ratios=ab_ratios,
        transcript_consequences=transcript_consequences,
    )


def create_structural_variant(var: 'cyvcf2.Variant', samples: list[str]):
    """
    takes an SV and creates a Model from it
    far less complicated than the SmallVariant model

    Args:
        var ():
        samples ():
    """

    # variant link for SVs is the RSID value
    info: dict[str, Any] = {x.lower(): y for x, y in var.INFO} | {'var_link': var.ID}

    # valid processing of inter-chromosomal SVs
    if all(attribute in info for attribute in ('chr2', 'end2')):
        info['svlen'] = f'{info["chr2"]}:{info["end2"]}'

    coordinates = Coordinates(chrom=var.CHROM.replace('chr', ''), pos=var.POS, ref=var.ALT[0], alt=str(info['svlen']))

    het_samples, hom_samples = get_non_ref_samples(variant=var, samples=samples)

    # set the class attributes
    boolean_categories = [key for key in info if key.startswith('categoryboolean')]

    # overwrite with true booleans
    for cat in boolean_categories:
        info[cat] = info.get(cat, 0) == 1

    return StructuralVariant(
        coordinates=coordinates,
        info=info,
        het_samples=het_samples,
        hom_samples=hom_samples,
        boolean_categories=boolean_categories,
    )


# CompHetDict structure: {sample: {variant_string: [variant, ...]}}
# sample: string, e,g, CGP12345
CompHetDict = dict[str, dict[str, list[VARIANT_MODELS]]]
GeneDict = dict[str, list[VARIANT_MODELS]]


def canonical_contigs_from_vcf(reader) -> set[str]:
    """
    read the header fields from the VCF handle
    return a set of all 'canonical' contigs

    Args:
        reader (cyvcf2.VCFReader):
    """

    # contig matching regex - remove all HLA/decoy/unknown
    contig_re = re.compile(r'^(chr)?[0-9XYMT]{1,2}$')

    return {
        contig['ID']
        for contig in reader.header_iter()
        if contig['HeaderType'] == 'CONTIG' and re.match(contig_re, contig['ID'])
    }


def gather_gene_dict_from_contig(
    contig: str,
    variant_source: 'cyvcf2.VCFReader',
    sv_source: 'cyvcf2.VCFReader | None' = None,
) -> GeneDict:
    """
    takes a cyvcf2.VCFReader instance, and a specified chromosome
    iterates over all variants in the region, and builds a lookup

    optionally takes a second VCF and incorporates into same dict

    Args:
        contig (): contig name from VCF header
        variant_source (): the VCF reader instance
        sv_source (): an optional list of SV VCFs

    Returns:
        A lookup in the form
        {
            gene1: [var1, var2],
            gene2: [var3],
            ...
        }
    """

    if bl_file := config_retrieve(['GeneratePanelData', 'blacklist'], ''):
        blacklist = read_json_from_path(bl_file, default=[])
    else:
        blacklist = []

    if not isinstance(blacklist, list):
        raise TypeError(f'Blacklist should be a list: {blacklist}')

    # a dict to allow lookup of variants on this whole chromosome
    contig_variants = 0
    contig_dict = defaultdict(list)

    # iterate over all variants on this contig and store by unique key
    # if contig has no variants, prints an error and returns []
    for variant in variant_source(contig):
        small_variant = create_small_variant(var=variant, samples=variant_source.samples)

        if small_variant is None:
            continue

        if small_variant.coordinates.string_format in blacklist:
            logger.info(f'Skipping blacklisted variant: {small_variant.coordinates.string_format}')
            continue

        # update the variant count
        contig_variants += 1

        # update the gene index dictionary
        contig_dict[small_variant.info.get('gene_id')].append(small_variant)

    # parse the SV VCF if provided, but not a necessary part of processing
    if sv_source:
        structural_variants = 0
        for variant in sv_source(contig):
            # create an abstract SV variant
            structural_variant = create_structural_variant(var=variant, samples=sv_source.samples)

            # update the variant count
            structural_variants += 1

            # update the gene index dictionary
            contig_dict[structural_variant.info.get('gene_id')].append(structural_variant)

        logger.info(f'Contig {contig} contained {structural_variants} SVs')

    logger.info(f'Contig {contig} contained {contig_variants} variants, in {len(contig_dict)} genes')

    return contig_dict


def read_json_from_path(read_path: str | None = None, default: Any = None, return_model: Any = None) -> Any:
    """
    take a path to a JSON file, read into an object
    if the path doesn't exist - return the default object
    uses cloudpath to be deployment agnostic

    Args:
        read_path (str): where to read from - if None... will return the value "default"
        default (Any):
        return_model (pydantic Models): any Model to read/validate as

    Returns:
        either the object from the JSON file, or None
    """

    if read_path is None:
        logger.error('read_json_from_path was passed the path "None"')
        return default

    assert isinstance(read_path, str)
    read_anypath = to_anypath(read_path)

    if not read_anypath.exists():
        logger.error(f'{read_path} did not exist')
        return default

    with read_anypath.open() as handle:
        json_data = json.load(handle)
        if return_model:
            # potentially walk-up model version
            model_data = lift_up_model_version(json_data, return_model)
            return return_model.model_validate(model_data)
        return json_data


def get_non_ref_samples(variant: 'cyvcf2.Variant', samples: list[str]) -> tuple[set[str], set[str]]:
    """
    For this variant, find all samples with a call. cyvcf2 uses 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT.
    Returns 2 sets of strings; het sample ID, hom sample IDs
    """
    het_samples = set()
    hom_samples = set()

    # this iteration is based on the cyvcf2 representations
    for sam, genotype_int in zip(samples, variant.gt_types, strict=True):
        if genotype_int in BAD_GENOTYPES:
            continue
        if genotype_int == HETALT:
            het_samples.add(sam)
        if genotype_int == HOMALT:
            hom_samples.add(sam)

    return het_samples, hom_samples


def extract_csq(csq_contents: str) -> list[dict]:
    """Handle extraction of the CSQ entries based on string in config."""

    # allow for no CSQ data, i.e. splice variant
    if not csq_contents:
        return []

    # break mono-CSQ-string into components
    csq_categories = config_retrieve(['RunHailFiltering', 'csq_string'])

    # iterate over all consequences, and make each into a dict
    txc_dict = [dict(zip(csq_categories, each_csq.split('|'), strict=True)) for each_csq in csq_contents.split(',')]

    # update this String to be either a float, or missing
    for each_dict in txc_dict:
        am_path = each_dict.get('am_pathogenicity')
        each_dict['am_pathogenicity'] = float(am_path) if am_path else ''

    return txc_dict


def find_comp_hets(var_list: list[VARIANT_MODELS], pedigree: PedigreeParser) -> CompHetDict:
    """
    Find compound het pairs, variants provided in the format [var1, var2, ...]

    generates pair content in the form
    {
        sample: {
            var_as_string: [partner_variant, ...],
        }
    }
    """

    # create an empty dictionary
    comp_het_results: CompHetDict = defaultdict(dict)

    # use combinations_with_replacement to find all gene pairs
    for var_1, var_2 in combinations_with_replacement(var_list, 2):
        assert var_1.coordinates.chrom == var_2.coordinates.chrom

        if (var_1.coordinates == var_2.coordinates) or var_1.coordinates.chrom in NON_HOM_CHROM:
            continue

        # iterate over any samples with a het overlap
        for sample in var_1.het_samples.intersection(var_2.het_samples):
            phased = False

            # don't assess male compound hets on sex chromosomes
            if pedigree.participants[sample].sex == 1 and var_1.coordinates.chrom in X_CHROMOSOME:
                continue

            # check for both variants being in the same phase set
            if sample in var_1.phased and sample in var_2.phased:
                # check for presence of the same phase set
                for phase_set in [ps for ps in var_1.phased[sample] if ps in var_2.phased[sample]]:
                    if var_1.phased[sample][phase_set] == var_2.phased[sample][phase_set]:
                        phased = True
            if not phased:
                comp_het_results[sample].setdefault(var_1.coordinates.string_format, []).append(var_2)
                comp_het_results[sample].setdefault(var_2.coordinates.string_format, []).append(var_1)

    return comp_het_results


def db_label_phenotypes(db_file: str, results: ResultData):
    """Using a SQLite DB, update the phenotype labels for all variants, annotate first phenotype match date."""
    connection = create_or_open_db(db_file)

    # get all current phenotype matches - this could be a series of lean queries, but we're expecting snowballing data
    # so we'd be querying for all samples anyway
    db_pheno_results = connection.execute(QUERY_ALL_VAR_STARS).fetchall()

    # index the results to be searchable
    pheno_dict = defaultdict(dict)
    variant_ids = {}
    for var_id, var_key, sample_id, first_pheno_match in db_pheno_results:
        pheno_dict[sample_id][var_key] = first_pheno_match
        variant_ids[var_key] = var_id

    new_var_stars = []
    update_var_stars = []

    for sample, var_list in results.results.items():
        sample_history = pheno_dict.get(sample, {})

        for var in var_list.variants:
            var_key = var.var_data.coordinates.string_format

            # pheno match date already exists in the db, must be an earlier date than today
            if early_date := sample_history.get(var_key):
                var.date_of_phenotype_match = early_date
                continue

            # no current phenotype match, nothing to add to the db
            if not var.phenotype_labels:
                continue

            # decide whether we need to insert or update
            if var_key in sample_history:
                # entry exists, but pheno match date was None (may have previously been ClinVar stars only) - update
                # this is the second round of history recording, no new variants should be added here
                update_var_stars.append((get_granular_date(), variant_ids[var_key], sample))
            else:
                # new entry, add it
                new_var_stars.append((variant_ids[var_key], sample, get_granular_date()))
    if update_var_stars:
        connection.executemany(UPDATE_VARSTARS_PHENO, update_var_stars)
    if new_var_stars:
        connection.executemany(INSERT_VARSTARS_PHENO, new_var_stars)
    if new_var_stars or update_var_stars:
        connection.commit()
    connection.close()


def date_from_string(filename: str) -> str:
    """
    Takes a filename, finds a date. Internally consistent with the way this codebase writes datetimes into file paths.
    """
    date_search = re.search(DATE_RE, filename)
    if date_search:
        return date_search.group()
    raise ValueError(f'No date found in {filename}')


def find_latest_file(results_folder: str, ext: str = 'json') -> str | None:
    """
    takes a directory of files, and finds the latest
    Args:
        results_folder (): local or remote folder, or don't call this method
        ext (): the type of files we're looking for

    Returns: most recent file path, or None
    """

    logger.info(f'Using results from {results_folder}')

    # this is currently a CloudPath to access globbing for files in cloud or local settings
    date_files = {}
    for filename in to_anypath(results_folder).glob(f'*.{ext}'):
        # if the filename is a valid date, add it to the dict, otherwise catch the error and skip it
        try:
            date_files[date_from_string(filename.name)] = filename
        except ValueError:
            logger.info(f'File {filename} did not have a valid date, skipping')

    if not date_files:
        logger.warning(f'The folder {results_folder} was provided, but did not contain any valid files')
        return None

    return str(date_files[max(date_files.keys())].absolute())


def generate_summary_stats(result_set: ResultData):
    """Reads over the variants present in the result set, and generates some per-sample and per-cohort stats."""

    # make this naive, i.e. not configured specifically based on a list of Categories in configuration
    category_count: dict[str, list[int]] = defaultdict(list)

    for sample_results in result_set.results.values():
        if len(sample_results.variants) == 0:
            continue

        sample_variants: dict[str, set[str]] = defaultdict(set)

        # iterate over all identified variants
        for each_var in sample_results.variants:
            var_string = each_var.var_data.coordinates.string_format

            # catch all comp-het pairs as a single variant - we are electing to count comp-het events as a single
            # variant, so we have to adjust our counting for this
            # do this by identifying all sorted pairwise combinations, sorting them, creating a single String for each,
            # and counting each unique String once
            if sups := each_var.support_vars:
                var_strings = [' '.join(sorted(combo)) for combo in combinations([var_string, *sups], 2)]

            # if the variant is dominant or Hom, count it once - list of length 1
            else:
                var_strings = [var_string]

            for unique_var in var_strings:
                # populate the 'any's
                sample_variants['any'].add(unique_var)

                # find all categories associated with this variant. For each category, add to corresponding list and set
                for category_value in each_var.categories:
                    sample_variants[category_value].add(unique_var)

        # record that these were the variants seen for this sample,
        for key, set_of_strings in sample_variants.items():
            category_count[key].append(len(set_of_strings))

    # update the counts-per-sample dictionary
    number_of_samples = len(result_set.results)

    stats_count: dict[str, dict[str, int | float]] = {}
    for category_type, count_list in category_count.items():
        # pad the observed counts to the number of samples under consideration
        count_list.extend([0] * (number_of_samples - len(count_list)))

        stats_count[category_type] = {
            'total': sum(count_list),
            'mean': sum(count_list) / number_of_samples,
            'min': min(count_list),
            'max': max(count_list),
            'median': statistics.median(count_list),
            'mode': statistics.mode(count_list) if len(count_list) > 0 else 0,
            'stddev': statistics.stdev(count_list) if len(count_list) > 1 else 0.0,
        }

    result_set.metadata.variant_breakdown = stats_count


def create_or_open_db(db_path: str) -> sqlite3.Connection:
    """
    Create or open a SQLite database at the specified path,

    Args:
        db_path (str): The file path for the SQLite database.

    Returns:
        sqlite3.Connection: A connection object to the SQLite database.
    """

    if to_anypath(db_path).exists():
        logger.info(f'Opening existing database at {db_path}')
        return sqlite3.connect(db_path)

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute('PRAGMA foreign_keys = ON;')

    # create a table to hold raw variant coord/ref/alt details
    variants_create = """
        CREATE TABLE IF NOT EXISTS variants (
            var_id INTEGER PRIMARY KEY,
            contig TEXT NOT NULL,
            position INTEGER NOT NULL,
            reference TEXT NOT NULL,
            alternate TEXT NOT NULL,
            UNIQUE (contig, position, reference, alternate)
        );
    """
    cursor.execute(variants_create)

    # record all variant categorisation events - sample, date, category
    categories_create = """
    CREATE TABLE IF NOT EXISTS categories (
        var_id INTEGER NOT NULL,
        sample_id TEXT NOT NULL,
        category TEXT NOT NULL,
        date TEXT NOT NULL,
        PRIMARY KEY (var_id, sample_id, category),
        FOREIGN KEY (var_id) REFERENCES variants(var_id),
        UNIQUE (var_id, sample_id, category)
    );
    """
    cursor.execute(categories_create)

    # record all variant stars and phenotype matches - dependent on sample & variant, but distinct from category
    var_star_create = """
    CREATE TABLE IF NOT EXISTS var_stars (
        var_id INTEGER NOT NULL,
        sample_id TEXT NOT NULL,
        clinvar_stars INTEGER,
        first_pheno_match TEXT,
        PRIMARY KEY (var_id, sample_id),
        FOREIGN KEY (var_id) REFERENCES variants(var_id),
        UNIQUE (var_id, sample_id)
    );
    """
    cursor.execute(var_star_create)

    # table to record comp-het partners
    partners_create = """
    CREATE TABLE IF NOT EXISTS partners (
        var_id INTEGER NOT NULL,
        sample_id TEXT NOT NULL,
        partners TEXT NOT NULL,
        PRIMARY KEY (var_id, sample_id),
        FOREIGN KEY (var_id) REFERENCES variants(var_id),
        UNIQUE (var_id, sample_id)
    );
    """
    cursor.execute(partners_create)
    conn.commit()

    return conn


def get_and_organise_db_results(sample_id: str, connection: sqlite3.Connection) -> dict:
    """ """
    history: dict[str, dict] = {}
    results = connection.cursor().execute(QUERY_SAMPLE_ALL.format(sample_id))
    for row in results.fetchall():
        (
            var_id,
            var_key,
            category,
            date,
            clinvar_stars,
            partners,
        ) = row

        # maybe we saw this before, with a different category attached
        if var_key in history:
            history[var_key]['categories'][category] = date
            continue

        # or maybe we didn't
        history[var_key] = {
            'var_id': var_id,
            'categories': {category: date},
            'clinvar_stars': clinvar_stars,
            'partners': partners.split(',') if partners else [],
        }
    return history


def get_all_db_variants(connection: sqlite3.Connection) -> dict[str, int]:
    """Query the DB for every variant we've seen before."""
    results = connection.cursor().execute(QUERY_ALL_VAR_IDS)
    return {var_key: var_id for var_id, var_key in results.fetchall()}


def db_date_annotate_results(current: ResultData, connection: sqlite3.Connection):
    """
    takes the current data, and annotates with previous dates if found
    build/update the historic data within the same loop.
    """

    # get all previously seen variants, and their IDs
    all_variant_ids = get_all_db_variants(connection=connection)

    new_categories = []
    new_partners = []
    new_var_stars = []
    update_var_stars = []

    # grab a cursor, treat yourself
    cursor = connection.cursor()

    # iterate over each sample section in the results
    for sample, content in current.results.items():
        # query the DB for results previously seen for this sample, can be none
        known_results = get_and_organise_db_results(sample, connection=connection)

        # check each variant found in this round
        for var in content.variants:
            # get the number of clinvar stars, if appropriate
            clinvar_stars = var.var_data.info.get('clinvar_stars', None)

            var_id = var.var_data.coordinates.string_format
            current_cats = var.categories

            # never seen this sample X variant before, everything is new
            if var_id not in known_results:
                # create the variant, store its ID
                cursor.execute(
                    INSERT_VAR_QUERY,
                    (
                        var.var_data.coordinates.chrom,
                        var.var_data.coordinates.pos,
                        var.var_data.coordinates.ref,
                        var.var_data.coordinates.alt,
                    ),
                )
                if last_row_id := cursor.lastrowid:
                    all_variant_ids[var_id] = last_row_id
                else:
                    connection.close()
                    raise ValueError(f'Failed to retrieve lastrowid after inserting new variant: {var_id}')

                for each_cat in current_cats:
                    new_categories.append(
                        (last_row_id, sample, translate_category(each_cat), get_granular_date()),
                    )

                if var.support_vars:
                    new_partners.append((last_row_id, sample, ','.join(sorted(var.support_vars))))

                if clinvar_stars is not None:
                    new_var_stars.append((last_row_id, sample, clinvar_stars))

                continue

            # we saw it before, so mark down the previously seen dates
            latest_date = sorted(known_results[var_id]['categories'].values(), reverse=True)[0]
            for each_cat in current_cats:
                if each_cat not in known_results[var_id]['categories']:
                    # add it to the list for creation
                    new_categories.append(
                        (
                            all_variant_ids[var_id],
                            sample,
                            translate_category(each_cat),
                            get_granular_date(),
                        ),
                    )

                # if we haven't seen this category before, update the latest date of evidence change
                else:
                    latest_date = get_granular_date()

            # update the number of stars if it's increased
            if clinvar_stars is not None:
                historic_stars = known_results[var_id]['clinvar_stars']
                if (historic_stars is None) or clinvar_stars > historic_stars:
                    var.clinvar_increase = True
                    update_var_stars.append((clinvar_stars, known_results[var_id]['var_id'], sample))

            # latest _new_ category date as evidence_last_changed timestamp
            var.evidence_last_updated = latest_date

    cursor.executemany(INSERT_CAT_QUERY, new_categories)
    cursor.executemany(INSERT_PARTNER_QUERY, new_partners)
    cursor.executemany(INSERT_VARSTARS_QUERY, new_var_stars)
    cursor.executemany(UPDATE_CLINVAR_QUERY, update_var_stars)
    connection.commit()
    connection.close()
