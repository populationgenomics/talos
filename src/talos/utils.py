"""
classes and methods shared across reanalysis components

HTTPX requests are backoff-wrapped using tenacity
https://tenacity.readthedocs.io/en/latest/
"""

import json
import re
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
from mendelbrot.bcftools_interpreter import TYPES_RE, classify_change
from mendelbrot.pedigree_parser import PedigreeParser
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
# don't fail, just suggest that someone raises an issue on github, but only print this message once
PHASE_BROKEN: bool = False

# not in use yet, if we want to reduce the number of TranscriptConsequence entries we store, this is a good start
BORING_CONSEQUENCES = ['downstream_gene_variant', 'intron_variant', 'upstream_gene_variant']


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
) -> SmallVariant | None:
    """Takes a small variant and creates a Model from it."""

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
    # keep it super flexible, lower case and a translated version
    ignored_categories = {word.lower() for word in config_retrieve(['ValidateMOI', 'ignore_categories'], [])}
    ignored_categories.update({translate_category(x) for x in ignored_categories})

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
    support_categories = {word.lower() for word in config_retrieve(['ValidateMOI', 'support_categories'], [])}
    support_categories.update({translate_category(x) for x in support_categories})

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

    # check that there's at least one category left after the PM5/Exomiser/SVDB/other processing, else return None
    # this isn't a sample-specific check, just a check that there's anything left worth classifying on
    if not (any(info[cat] for cat in boolean_categories) or any(info[cat] for cat in sample_categories)):
        return None

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


def create_structural_variant(var: 'cyvcf2.Variant', samples: list[str]) -> StructuralVariant:
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
    contig_re = re.compile(r'^(chr)?[0-9XY]{1,2}$')

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
        if (small_variant := create_small_variant(var=variant, samples=variant_source.samples)) is None:
            continue

        if small_variant.coordinates.string_format in blacklist:
            logger.info(f'Skipping blacklisted variant: {small_variant.coordinates.string_format}')
            continue

        # update the variant count
        contig_variants += 1

        # update the gene index dictionary
        contig_dict[small_variant.info['gene_id']].append(small_variant)

    # parse the SV VCF if provided, but not a necessary part of processing
    if sv_source:
        structural_variants = 0
        for variant in sv_source(contig):
            # create an abstract SV variant
            structural_variant = create_structural_variant(var=variant, samples=sv_source.samples)

            # update the variant count
            structural_variants += 1

            # update the gene index dictionary
            contig_dict[structural_variant.info['gene_id']].append(structural_variant)

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
    txc_dicts = [dict(zip(csq_categories, each_csq.split('|'), strict=True)) for each_csq in csq_contents.split(',')]

    # update this String to be either a float, or missing
    for each_dict in txc_dicts:
        am_path = each_dict.get('am_pathogenicity')
        each_dict['am_pathogenicity'] = float(am_path) if am_path else ''

        # update the BCFtools amino acid change with an HGVS reinterpretation
        csq_string = str(each_dict['consequence'])
        if (type_match := TYPES_RE.match(csq_string)) and (aa := each_dict.get('amino_acid_change')):
            each_dict['amino_acid_change'] = classify_change(
                aa,
                consequence=type_match[0],
            )

    return txc_dicts


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


def annotate_variant_dates_using_prior_results(results: ResultData, previous_results: ResultData | None = None):
    """
    loads the most recent prior result set (if it exists)
    annotates previously seen variants with the most recent date seen
    write two files (total, and latest - previous)

    Args:
        results (ResultData): the results produced during this run
        previous_results (ResultData | None): optional, results of a previous analysis to build on

    Returns:
        None, but if previous results are provided, updates the results in place

    """

    if previous_results is None:
        logger.info('No previous results provided, saving with current dates only')
        return

    # for every sample-variant in both datasets, check for new categories, and adjust all dates accordingly
    # bonus behaviour, if a prev. variant is not present in the current results, insert with found_in_current_run=False
    for sample, content in results.results.items():
        # this variant was absent from previous analysis, or was newly added in this run
        if sample not in previous_results.results:
            continue

        # collect variant list into a dict for easy lookup
        new_vars = {var.var_data.coordinates.string_format: var for var in content.variants}

        # scroll over all the variants found in the previous run
        for old_var in previous_results.results[sample].variants:
            old_coord = old_var.var_data.coordinates.string_format
            # we didn't recover this variant in the current run - insert it into new results, marking as not found
            if old_coord not in new_vars:
                old_var.categories = {translate_category(key): val for key, val in old_var.categories.items()}
                old_var.found_in_current_run = False
                old_var.clinvar_increase = False
                content.variants.append(old_var)
                continue

            # we found this variant again, but do we have any new categories to add?
            new_var = new_vars[old_coord]
            if new_var.clinvar_stars:
                new_var.clinvar_increase = bool(
                    old_var.clinvar_stars is None or new_var.clinvar_stars > old_var.clinvar_stars,
                )

            for cat, date in old_var.categories.items():
                new_var.categories[translate_category(cat)] = date

            # collect all the dates we have for first category assignment
            category_dates = list(new_var.categories.values())

            # we previously had a phenotype match date, carry it forward
            if old_pheno := old_var.date_of_phenotype_match:
                new_var.date_of_phenotype_match = old_pheno
                category_dates.append(old_pheno)

            # new supporting comp-het partners = new evidence change date
            if new_var.support_vars - old_var.support_vars:
                category_dates.append(get_granular_date())
            new_var.support_vars.update(old_var.support_vars)

            new_var.evidence_last_updated = max(category_dates)
            new_var.first_tagged = min(new_var.categories.values())


def generate_summary_stats(result_set: ResultData):
    """Reads over the variants present in the result set, and generates some per-sample and per-cohort stats."""

    # make this naive, i.e. not configured specifically based on a list of Categories in configuration
    category_count: dict[str, list[int]] = defaultdict(list)
    samples_with_no_variants: list[str] = []

    for sample_id, sample_results in result_set.results.items():
        if len(sample_results.variants) == 0:
            samples_with_no_variants.append(sample_id)
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

    # record samples with no variants (by sample ID)
    result_set.metadata.samples_with_no_variants = samples_with_no_variants

    # Compute unused external labels, if an external labels file is configured
    unused_ext_labels: list[dict[str, Any]] = []
    labels_path = config_retrieve(['CreateTalosHTML', 'external_labels'], None)
    if labels_path:
        ext_label_map: dict[str, dict[str, list[str]]] = read_json_from_path(labels_path, default={}) or {}

        # Remove any labels that correspond to observed variants
        for sample_id, sample_results in result_set.results.items():
            if sample_id not in ext_label_map:
                continue
            for variant in sample_results.variants:
                var_string = variant.var_data.coordinates.string_format
                ext_label_map[sample_id].pop(var_string, None)

        # Whatever remains in ext_label_map did not match any observed variants
        for sam, var_dict in ext_label_map.items():
            for var_id, labels in var_dict.items():
                unused_ext_labels.append({'sample': sam, 'variant': var_id, 'labels': labels})

    result_set.metadata.unused_ext_labels = unused_ext_labels
