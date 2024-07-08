"""
classes and methods shared across reanalysis components
"""

import json
import re
from collections import defaultdict
from datetime import datetime
from itertools import chain, combinations_with_replacement, islice
from pathlib import Path
from string import punctuation
from typing import Any

import backoff
import cyvcf2
import requests
import zoneinfo
from backoff import fibo
from cloudpathlib.anypath import to_anypath
from peds import open_ped
from requests.exceptions import ReadTimeout, RequestException

from talos.config import config_retrieve
from talos.models import (
    VARIANT_MODELS,
    CategoryMeta,
    Coordinates,
    FileTypes,
    HistoricPanels,
    HistoricSampleVariant,
    HistoricVariants,
    PanelApp,
    Pedigree,
    PedigreeMember,
    PhenoPacketHpo,
    PhenotypeMatchedPanels,
    ResultData,
    SmallVariant,
    StructuralVariant,
    lift_up_model_version,
)
from talos.static_values import get_granular_date, get_logger

HOMREF: int = 0
HETALT: int = 1
UNKNOWN: int = 2
HOMALT: int = 3

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
REMOVE_IN_SINGLETONS = {'categorysample4'}

DATE_RE = re.compile(r'\d{4}-\d{2}-\d{2}')


def make_flexible_pedigree(pedigree: str) -> Pedigree:
    """
    takes the representation offered by peds and reshapes it to be searchable
    this is really just one short step from writing my own implementation...

    Args:
        pedigree (str): path to a pedigree file

    Returns:
        a searchable representation of the ped file
    """
    new_ped = Pedigree()
    ped_data = open_ped(pedigree)
    for family in ped_data:
        for member in family:
            # any extra columns are parsed as a tuple of strings
            member_data = member.data

            # default to repeating internal ID
            ext_id = member_data[0] if member_data else member.id

            # can be an empty list
            hpos = [PhenoPacketHpo(id=hpo, label=hpo) for hpo in member_data[1:]]

            me = PedigreeMember(
                family=member.family,
                id=member.id,
                mother=None if member.mom == '0' else member.mom,
                father=None if member.dad == '0' else member.dad,
                sex=member.sex,
                affected=member.phenotype,
                ext_id=ext_id,
                hpo_terms=hpos,
            )

            # add as a member
            new_ped.members.append(me)

            # add to a lookup
            new_ped.by_id[me.id] = me

            # add to a list of members in this family
            new_ped.by_family.setdefault(me.family, []).append(me)
    return new_ped


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


def identify_file_type(file_path: str) -> FileTypes | Exception:
    """
    return type of the file, if present in FileTypes enum

    Args:
        file_path (str):

    Returns:
        A matching file type, or die
    """
    pl_filepath = Path(file_path)

    # pull all extensions (e.g. .vcf.bgz will be split into [.vcf, .bgz]
    extensions = pl_filepath.suffixes

    assert len(extensions) > 0, 'cannot identify input type from extensions'

    if extensions[-1] == '.ht':
        return FileTypes.HAIL_TABLE
    if extensions[-1] == '.mt':
        return FileTypes.MATRIX_TABLE
    if extensions == ['.vcf']:
        return FileTypes.VCF
    if extensions == ['.vcf', '.gz']:
        return FileTypes.VCF_GZ
    if extensions == ['.vcf', '.bgz']:
        return FileTypes.VCF_BGZ
    raise TypeError(f'File cannot be definitively typed: {str(extensions)}')


@backoff.on_exception(
    wait_gen=fibo,
    exception=(TimeoutError, ReadTimeout, RequestException),
    max_time=200,
    logger=get_logger(),
)
def get_json_response(url):
    """
    takes a request URL, checks for healthy response, returns the JSON
    For this purpose we only expect a dictionary return
    List use-case (activities endpoint) no longer supported

    Args:
        url (str): URL to retrieve JSON format data from

    Returns:
        the JSON response from the endpoint
    """

    response = requests.get(url, headers={'Accept': 'application/json'}, timeout=60)
    response.raise_for_status()  # Raise an exception for bad responses (4xx and 5xx)
    return response.json()


def get_new_gene_map(
    panelapp_data: PanelApp,
    pheno_panels: PhenotypeMatchedPanels | None = None,
) -> dict[str, set[str]]:
    """
    The aim here is to generate a list of all the samples for whom
    a given gene should be treated as new during this analysis. This
    prevents the need for back-filtering results at the end of
    classification.

    Generate a map of
    {gene: [samples, where, this, is, 'new']}

    Args:
        panelapp_data ():
        pheno_panels (PhenotypeMatchedPanels):
    """

    # any dataset-specific panel data, + 'core' panel
    cohort_panels = [
        *config_retrieve(['GeneratePanelData', 'forced_panels'], []),
        config_retrieve(['GeneratePanelData', 'default_panel']),
    ]

    # collect all genes new in at least one panel
    new_genes: dict[str, set[int]] = {ensg: content.new for ensg, content in panelapp_data.genes.items() if content.new}

    # if there's no panel matching, new applies to everyone
    if pheno_panels is None:
        return {ensg: {'all'} for ensg in new_genes}

    # if we have pheno-matched participants, more complex
    panel_samples: dict[int, set[str]] = defaultdict(set)

    # double layered iteration, but only on a small object
    for sample, data in pheno_panels.samples.items():
        for panel in data.panels:
            panel_samples[panel].add(sample)

    pheno_matched_new = {}

    # iterate over the new genes and find out who they are new for
    for gene, panels in new_genes.items():
        if any(panel in cohort_panels for panel in panels):
            pheno_matched_new[gene] = {'all'}
            continue

        # else, find the specific samples
        samples = set()
        for panel_id in panels:
            # this line causes problems running mismatched pheno/panel data
            if panel_id not in panel_samples:
                raise AssertionError(f'PanelID {panel_id} not attached to any samples')
            samples.update(panel_samples[panel_id])
        pheno_matched_new[gene] = samples

    return pheno_matched_new


def get_phase_data(samples, var) -> dict[str, dict[int, str]]:
    """
    read phase data from this variant

    Args:
        samples ():
        var ():
    """
    phased_dict: dict[str, dict[int, str]] = defaultdict(dict)

    # first set the numpy.ndarray to be a list of ints
    # then zip against ordered sample IDs
    # this might need to store the exact genotype too
    # i.e. 0|1 and 1|0 can be in the same phase-set
    # but are un-phased variants

    try:
        for sample, phase, genotype in zip(samples, map(int, var.format('PS')), var.genotypes):
            # cyvcf2.Variant holds two ints, and a bool
            allele_1, allele_2, phased = genotype
            if not phased:
                continue
            gt = f'{allele_1}|{allele_2}'
            # phase set is a number
            if phase != PHASE_SET_DEFAULT:
                phased_dict[sample][phase] = gt
    except KeyError as ke:
        get_logger().info('failed to find PS phase attributes')
        try:
            # retry using PGT & PID
            for sample, phase_gt, phase_id in zip(samples, var.format('PGT'), var.format('PID')):
                if phase_gt != '.' and phase_id != '.':
                    phased_dict[sample][phase_id] = phase_gt

        except KeyError as ke2:
            get_logger().info('also failed using PID and PGT')
            raise ke from ke2

    return dict(phased_dict)


def organise_pm5(info_dict: dict[str, Any]) -> dict[str, Any]:
    """
    method dedicated to handling the new pm5 annotations

    e.g. categorydetailsPM5=27037::Pathogenic::1+27048::Pathogenic::1;
    1. break into component allele data

    Returns:
        None, updates self. attributes
    """

    if 'categorydetailspm5' not in info_dict:
        return info_dict

    pm5_content = info_dict.pop('categorydetailspm5')

    # nothing to do here
    if pm5_content == 'missing':
        info_dict['categorybooleanpm5'] = 0
        return info_dict

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

    return info_dict


def create_small_variant(
    var: cyvcf2.Variant,
    samples: list[str],
    as_singletons=False,
    new_genes: dict[str, set[str]] | None = None,
):
    """
    takes a small variant and creates a Model from it

    Args:
        var ():
        samples ():
        as_singletons ():
        new_genes ():
    """

    coordinates = Coordinates(chrom=var.CHROM.replace('chr', ''), pos=var.POS, ref=var.REF, alt=var.ALT[0])
    depths: dict[str, int] = dict(zip(samples, map(int, var.gt_depths)))
    info: dict[str, Any] = {x.lower(): y for x, y in var.INFO} | {'seqr_link': coordinates.string_format}

    # optionally - ignore some categories from this analysis
    if ignore_cats := config_retrieve(['ValidateMOI', 'ignore_categories'], []):
        info = {key: val for key, val in info.items() if key not in ignore_cats}

    het_samples, hom_samples = get_non_ref_samples(variant=var, samples=samples)

    # hot-swap cat 2 from a boolean to a sample list - if appropriate
    if info.get('categoryboolean2', 0) and new_genes:
        new_gene_samples: set[str] = new_genes.get(info.get('gene_id', 'MISSING'), set())

        # if 'all', keep cohort-wide boolean flag
        if new_gene_samples == {'all'}:
            get_logger().debug('New applies to all samples')

        # otherwise assign only a specific sample list
        elif new_gene_samples:
            _boolcat = info.pop('categoryboolean2')
            info['categorysample2'] = new_gene_samples

    # organise PM5
    info = organise_pm5(info)

    # set the class attributes
    boolean_categories = [key for key in info if key.startswith('categoryboolean')]
    sample_categories = [key for key in info if key.startswith('categorysample')]
    support_categories = [key for key in info if key.startswith('categorysupport')]

    # overwrite with true booleans
    for cat in support_categories + boolean_categories:
        info[cat] = info.get(cat, 0) == 1

    # sample categories are a set of strings or 'missing'
    # if cohort runs as singletons, remove possibility of de novo
    # if not singletons, split each into a set of sample IDs
    for sam_cat in sample_categories:
        if as_singletons and sam_cat in REMOVE_IN_SINGLETONS:
            info[sam_cat] = set()
        elif isinstance(info[sam_cat], str):
            info[sam_cat] = info[sam_cat].split(',') if info[sam_cat] != 'missing' else set()
        elif isinstance(info[sam_cat], list):
            info[sam_cat] = set(info[sam_cat])

    phased = get_phase_data(samples, var)
    ab_ratios = dict(zip(samples, map(float, var.gt_alt_freqs)))
    transcript_consequences = extract_csq(csq_contents=info.pop('csq', ''))

    return SmallVariant(
        coordinates=coordinates,
        info=info,
        het_samples=het_samples,
        hom_samples=hom_samples,
        boolean_categories=boolean_categories,
        sample_categories=sample_categories,
        sample_support=support_categories,
        phased=phased,
        depths=depths,
        ab_ratios=ab_ratios,
        transcript_consequences=transcript_consequences,
    )


def create_structural_variant(var: cyvcf2.Variant, samples: list[str]):
    """
    takes an SV and creates a Model from it
    far less complicated than the SmallVariant model

    Args:
        var ():
        samples ():
    """

    info: dict[str, Any] = {x.lower(): y for x, y in var.INFO}

    # valid processing of inter-chromosomal SVs
    if all(attribute in info for attribute in ('chr2', 'end2')):
        info['svlen'] = f'{info["chr2"]}:{info["end2"]}'

    # this is the right ID for Seqr
    info['seqr_link'] = info['variantid']

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
    variant_source,
    new_gene_map: dict[str, set[str]],
    singletons: bool = False,
    sv_sources: list | None = None,
) -> GeneDict:
    """
    takes a cyvcf2.VCFReader instance, and a specified chromosome
    iterates over all variants in the region, and builds a lookup

    optionally takes a second VCF and incorporates into same dict

    Args:
        contig (): contig name from VCF header
        variant_source (): the VCF reader instance
        sv_sources (): an optional list of SV VCFs
        new_gene_map ():
        singletons ():

    Returns:
        A lookup in the form
        {
            gene1: [var1, var2],
            gene2: [var3],
            ...
        }
    """
    if sv_sources is None:
        sv_sources = []
    if bl_file := config_retrieve(['GeneratePanelData', 'blacklist'], ''):
        blacklist = read_json_from_path(bl_file, default=[])
    else:
        blacklist = []

    assert isinstance(blacklist, list)

    # a dict to allow lookup of variants on this whole chromosome
    contig_variants = 0
    contig_dict = defaultdict(list)

    # iterate over all variants on this contig and store by unique key
    # if contig has no variants, prints an error and returns []
    for variant in variant_source(contig):
        small_variant = create_small_variant(
            var=variant,
            samples=variant_source.samples,
            as_singletons=singletons,
            new_genes=new_gene_map,
        )

        if small_variant.coordinates.string_format in blacklist:
            get_logger().info(f'Skipping blacklisted variant: {small_variant.coordinates.string_format}')
            continue

        # if unclassified, skip the whole variant
        if not small_variant.is_classified:
            continue

        # update the variant count
        contig_variants += 1

        # update the gene index dictionary
        contig_dict[small_variant.info.get('gene_id')].append(small_variant)

    for sv_source in sv_sources:
        structural_variants = 0
        for variant in sv_source(contig):
            # create an abstract SV variant
            structural_variant = create_structural_variant(var=variant, samples=sv_source.samples)

            # update the variant count
            structural_variants += 1

            # update the gene index dictionary
            contig_dict[structural_variant.info.get('gene_id')].append(structural_variant)

        get_logger().info(f'Contig {contig} contained {structural_variants} SVs')

    get_logger().info(f'Contig {contig} contained {contig_variants} variants, in {len(contig_dict)} genes')

    return contig_dict


def read_json_from_path(read_path: str | None = None, default: Any = None, return_model: Any = None) -> Any:
    """
    take a path to a JSON file, read into an object
    if the path doesn't exist - return the default object
    uses cloudpath to be deployment agnostic

    Args:
        read_path (str): where to read from - if None... will return default
        default (Any):
        return_model (pydantic Models): any Model to read/validate as

    Returns:
        either the object from the JSON file, or None
    """

    if read_path is None:
        get_logger().error('read_json_from_path was passed the path "None"')
        return default

    assert isinstance(read_path, str)
    read_anypath = to_anypath(read_path)

    if not read_anypath.exists():
        get_logger().error(f'{read_path} did not exist')
        return default

    with read_anypath.open() as handle:
        json_data = json.load(handle)
        if return_model:
            # potentially walk-up model version
            model_data = lift_up_model_version(json_data, return_model)
            return return_model.model_validate(model_data)
        return json_data


def get_simple_moi(input_mois: set[str], chrom: str) -> set[str]:
    """
    takes the vast range of PanelApp MOIs, and reduces to a
    range of cases which can be easily implemented in RD analysis

    Args:
        input_mois (set[str]): all the MOIs for this gene
        chrom ():
    """

    default = 'Hemi_Bi_In_Female' if chrom in X_CHROMOSOME else 'Biallelic'

    return_mois: set[str] = set()

    for input_moi in input_mois:
        # skip over ignore-able MOIs
        if input_moi in IRRELEVANT_MOI:
            continue

        # split each PanelApp MOI into a list of strings
        input_list = input_moi.translate(str.maketrans('', '', punctuation)).split()

        # run a match: case to classify it
        match input_list:
            case ['biallelic', *_additional]:
                return_mois.add('Biallelic')
            case ['both', *_additional]:
                return_mois.add('Mono_And_Biallelic')
            case ['monoallelic', *_additional]:
                if chrom in X_CHROMOSOME:
                    return_mois.add('Hemi_Mono_In_Female')
                else:
                    return_mois.add('Monoallelic')
            case ['xlinked', *additional] if 'biallelic' in additional:
                return_mois.add('Hemi_Bi_In_Female')
            case ['xlinked', *_additional]:
                return_mois.add('Hemi_Mono_In_Female')
            case _:
                continue

    # adda default - solves the all-irrelevant or empty-input cases
    if not return_mois:
        return_mois.add(default)

    return return_mois


def get_non_ref_samples(variant, samples: list[str]) -> tuple[set[str], set[str]]:
    """
    for this variant, find all samples with a call
    cyvcf2 uses 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT

    Args:
        variant (cyvcf2.Variant):
        samples (list[str]):

    Returns:
        2 sets of strings; het and hom
    """
    het_samples = set()
    hom_samples = set()

    # this iteration is based on the cyvcf2 representations
    for sam, genotype_int in zip(samples, variant.gt_types):
        if genotype_int in BAD_GENOTYPES:
            continue
        if genotype_int == HETALT:
            het_samples.add(sam)
        if genotype_int == HOMALT:
            hom_samples.add(sam)

    return het_samples, hom_samples


def extract_csq(csq_contents: str) -> list[dict]:
    """
    handle extraction of the CSQ entries based on string in config

    Args:
        csq_contents ():
    """

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


def find_comp_hets(var_list: list[VARIANT_MODELS], pedigree: Pedigree) -> CompHetDict:
    """
    manual implementation to find compound hets
    variants provided in the format

    [var1, var2, ..]

    generate pair content in the form
    {sample: {var_as_string: [partner_variant, ...]}}

    Args:
        var_list (list[VARIANT_MODELS]): all variants in this gene
        pedigree (): Pedigree
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
            if pedigree.by_id[sample].sex == '1' and var_1.coordinates.chrom in X_CHROMOSOME:
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


def generate_fresh_latest_results(current_results: ResultData, prefix: str = ''):
    """
    This will be called if a cohort has no latest results, but has
    indicated a place to save them.
    Args:
        current_results (ResultData): results from this current run
        prefix (str): optional prefix for the filename
    """

    new_history = HistoricVariants(metadata=CategoryMeta(categories=config_retrieve('categories', {})))
    for sample, content in current_results.results.items():
        for var in content.variants:
            # bank the number of clinvar stars, if any
            if '1' in var.categories:
                clinvar_stars = var.var_data.info.get('clinvar_stars')
                assert isinstance(clinvar_stars, int)
            else:
                clinvar_stars = None

            new_history.results.setdefault(sample, {})[var.var_data.coordinates.string_format] = HistoricSampleVariant(
                categories={cat: var.first_tagged for cat in var.categories},
                support_vars=var.support_vars,
                independent=var.independent,
                first_tagged=var.first_tagged,  # this could be min of available values, but this is first analysis
                clinvar_stars=clinvar_stars,
            )
    save_new_historic(results=new_history, prefix=prefix)


def filter_results(results: ResultData, singletons: bool):
    """
    loads the most recent prior result set (if it exists)
    annotates previously seen variants with the most recent date seen
    write two files (total, and latest - previous)

    Args:
        results (ResultData): the results produced during this run
        singletons (bool): whether to read/write a singleton specific file

    Returns: same results annotated with date-first-seen
    """

    historic_folder = config_retrieve('result_history')

    if historic_folder is None:
        get_logger().info('No historic data folder, no filtering')
        # update all the evidence_last_updated
        for content in results.results.values():
            for var in content.variants:
                var.evidence_last_updated = get_granular_date()
        return

    # get the latest result file from the folder
    # this will be none if the folder doesn't exist or is empty
    prefix = 'singletons_' if singletons else ''

    # 2 is the required prefix, i.e. 2022_*, to discriminate vs. 'singletons_'
    latest_results_path = find_latest_file(results_folder=config_retrieve('result_history', None), start=prefix or '2')

    get_logger().info(f'latest results: {latest_results_path}')

    # get latest results as a HistoricVariants object, or fail - on fail, return
    if (latest_results := read_json_from_path(latest_results_path, return_model=HistoricVariants)) is None:
        # generate and write some new latest data
        generate_fresh_latest_results(current_results=results, prefix=prefix)
        return

    assert isinstance(latest_results, HistoricVariants)

    date_annotate_results(results, latest_results)
    save_new_historic(results=latest_results, prefix=prefix)


def save_new_historic(results: HistoricVariants | HistoricPanels, prefix: str = ''):
    """
    save the new results in the historic results dir

    Args:
        results (HistoricVariants): object to save as JSON
        prefix (str): name prefix for this file (optional)
    """

    if (directory := config_retrieve('result_history', None)) is None:
        get_logger().info('No historic data folder, no saving')
        return

    # we're using cloud paths here
    new_file = to_anypath(directory) / f'{prefix}{TODAY}.json'
    with new_file.open('w') as handle:
        handle.write(results.model_dump_json(indent=4))

    get_logger().info(f'Wrote new data to {new_file}')


def date_from_string(string: str) -> str:
    """
    takes a string, finds the date. Simples
    Args:
        string (a filename):

    Returns:
        the String YYYY-MM-DD
    """
    date_search = re.search(DATE_RE, string)
    if date_search:
        return date_search.group()
    raise ValueError(f'No date found in {string}')


def find_latest_file(results_folder: str, start: str = '', ext: str = 'json') -> str | None:
    """
    takes a directory of files, and finds the latest
    Args:
        results_folder (): local or remote folder, or don't call this method
        start (str): the start of the filename, if applicable
        ext (): the type of files we're looking for

    Returns: most recent file path, or None
    """

    get_logger().info(f'Using results from {results_folder}')

    # this is currently a CloudPath to access globbing for files in cloud or local settings
    date_files = {
        date_from_string(filename.name): filename for filename in to_anypath(results_folder).glob(f'{start}*.{ext}')
    }
    if not date_files:
        return None

    return str(date_files[max(date_files.keys())].absolute())


def date_annotate_results(current: ResultData, historic: HistoricVariants):
    """
    takes the current data, and annotates with previous dates if found
    build/update the historic data within the same loop

    logically we can't reach this method in production without historic
    data being present. This method no longer permits historic data to be
    missing, and edits objects in place

    first_seen date is the earliest date of any category assignment
    this date should be static

    Args:
        current (ResultData): results generated during this run
        historic (HistoricVariants): optionally, historic data

    Returns:
        updated/instantiated cumulative data
    """

    # update to latest category descriptions
    historic.metadata.categories.update(config_retrieve('categories'))

    for sample, content in current.results.items():
        # get the historical record for this sample
        sample_historic = historic.results.get(sample, {})

        # check each variant found in this round
        for var in content.variants:
            # get the number of clinvar stars, if appropriate
            if '1' in var.categories:
                clinvar_stars = var.var_data.info.get('clinvar_stars')
                assert isinstance(clinvar_stars, int)
            else:
                clinvar_stars = None

            var_id = var.var_data.coordinates.string_format
            current_cats = var.categories

            # this variant was previously seen
            if hist := sample_historic.get(var_id):
                # bool if this was ever independent (i.e. dominant or Homozygous)
                # changed logic - do not update dates of prior category assignments
                if var.independent and not hist.independent:
                    hist.independent = True

                # if we have any new categories don't alter the date
                if new_cats := current_cats - set(hist.categories):
                    # add any new categories
                    for cat in new_cats:
                        hist.categories[cat] = get_granular_date()

                # update to include all support_variants
                hist.support_vars.update(var.support_vars)

                # mark the first seen timestamp
                var.first_tagged = hist.first_tagged

                # log an increase in ClinVar star rating
                if clinvar_stars:
                    historic_stars = hist.clinvar_stars

                    # if clinvar_stars was previously None, or was previously lower, store and keep a boolean
                    # then update the history for this variant to flag that the star rating has increased
                    if (historic_stars is None) or clinvar_stars > historic_stars:
                        var.clinvar_increase = True
                        hist.clinvar_stars = clinvar_stars

                # latest _new_ category date as evidence_last_changed timestamp
                var.evidence_last_updated = sorted(hist.categories.values(), reverse=True)[0]

            # totally new variant
            else:
                historic.results.setdefault(sample, {})[var_id] = HistoricSampleVariant(
                    categories={cat: get_granular_date() for cat in current_cats},
                    support_vars=var.support_vars,
                    independent=var.independent,
                    first_tagged=get_granular_date(),
                    clinvar_stars=clinvar_stars,
                )
