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
import peddy
import requests
from backoff import expo
from gql.gql import DocumentNode
from gql.transport.exceptions import TransportQueryError, TransportServerError

from cpg_utils import Path as CPGPathType
from cpg_utils import to_path
from cpg_utils.config import get_config
from metamist.graphql import query

from reanalysis.models import (
    VARIANT_MODELS,
    CategoryMeta,
    Coordinates,
    FileTypes,
    HistoricPanels,
    HistoricSampleVariant,
    HistoricVariants,
    PanelApp,
    PhenotypeMatchedPanels,
    ResultData,
    SmallVariant,
    StructuralVariant,
)
from reanalysis.static_values import get_granular_date, get_logger

HOMREF: int = 0
HETALT: int = 1
UNKNOWN: int = 2
HOMALT: int = 3

# in cyVCF2, these ints represent HOMREF, and UNKNOWN
BAD_GENOTYPES: set[int] = {HOMREF, UNKNOWN}
PHASE_SET_DEFAULT = -2147483648
NON_HOM_CHROM = ['X', 'Y', 'MT', 'M']
CHROM_ORDER = list(map(str, range(1, 23))) + NON_HOM_CHROM
X_CHROMOSOME = {'X'}
TODAY = datetime.now().strftime('%Y-%m-%d_%H:%M')

# most lenient to most conservative
# usage = if we have two MOIs for the same gene, take the broadest
ORDERED_MOIS = [
    'Mono_And_Biallelic',
    'Monoallelic',
    'Hemi_Mono_In_Female',
    'Hemi_Bi_In_Female',
    'Biallelic',
]
IRRELEVANT_MOI = {'unknown', 'other'}
REMOVE_IN_SINGLETONS = {'categorysample4'}

# global config holders
COHORT_CONFIG: dict | None = None
COHORT_SEQ_CONFIG: dict | None = None


@backoff.on_exception(wait_gen=expo, exception=(TransportQueryError, TransportServerError), max_time=20)
def wrapped_gql_query(query_node: DocumentNode, variables: dict[str, Any] | None = None) -> dict[str, Any]:
    """
    wrapped gql query method, with retries
    uses an exponential backoff retry timer to space out attempts

    Args:
        query_node (the result of a gql() call):
        variables (dict of parameters, or None):

    Returns:
        the response from the query
    """
    return query(query_node, variables=variables)


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


@backoff.on_exception(wait_gen=expo, exception=(requests.RequestException, TimeoutError), max_time=20)
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


def get_cohort_config(dataset: str | None = None):
    """
    return the cohort-specific portion of the config file, or fail

    Returns:
        the dict of cohort and genome/exome specific content
    """

    global COHORT_CONFIG
    if COHORT_CONFIG is None:
        dataset = dataset or get_config()['workflow']['dataset']
        COHORT_CONFIG = get_config().get('cohorts', {}).get(dataset)
        if COHORT_CONFIG is None:
            raise AssertionError(f'{dataset} is not represented in config')
    return COHORT_CONFIG


def get_cohort_seq_type_conf(dataset: str | None = None):
    """
    return the cohort-specific portion of the config file,
    chased down to the exome/genome specific portion

    Args:
        dataset (str): the dataset to retrieve config for
                       defaults to config/workflow/dataset

    Returns:
        the dict of cohort and genome/exome specific content
    """
    global COHORT_SEQ_CONFIG
    if COHORT_SEQ_CONFIG is None:
        dataset = dataset or get_config()['workflow']['dataset']
        cohort_conf = get_cohort_config(dataset)
        seq_type = get_config()['workflow']['sequencing_type']
        COHORT_SEQ_CONFIG = cohort_conf.get(seq_type, {})
        assert COHORT_SEQ_CONFIG, f'{dataset} - {seq_type} is not represented in config'
    return COHORT_SEQ_CONFIG


# todo rethink this whole method?
def get_new_gene_map(
    panelapp_data: PanelApp,
    pheno_panels: PhenotypeMatchedPanels | None = None,
    dataset: str | None = None,
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
        dataset ():

    Returns:

    """

    # any dataset-specific panel data, + 'core' panel
    cohort_panels = [
        *get_cohort_config(dataset).get('cohort_panels', []),
        get_config()['panels']['default_panel'],
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
        except KeyError:
            get_logger().info('also failed using PID and PGT')
            raise ke

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
    if ignore_cats := get_config()['workflow'].get('ignore_categories'):
        info = {key: val for key, val in info.items() if key not in ignore_cats}

    het_samples, hom_samples = get_non_ref_samples(variant=var, samples=samples)

    # hot-swap cat 2 from a boolean to a sample list - if appropriate
    if info.get('categoryboolean2', 0) and new_genes:
        new_gene_samples: set[str] = new_genes.get(info.get('gene_id'), set())  # type: ignore

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
    boolean_categories = [key for key in info.keys() if key.startswith('categoryboolean')]
    sample_categories = [key for key in info.keys() if key.startswith('categorysample')]
    support_categories = [key for key in info.keys() if key.startswith('categorysupport')]

    # overwrite with true booleans
    for cat in support_categories + boolean_categories:
        info[cat] = info.get(cat, 0) == 1

    # sample categories are a set of strings or 'missing'
    # if cohort runs as singletons, remove possibility of de novo
    # if not singletons, split each into a set of sample IDs
    # todo I have messed with this, check it works
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

    coordinates = Coordinates(
        chrom=var.CHROM.replace('chr', ''),
        pos=var.POS,
        ref=var.ALT[0],
        alt=str(info['svlen']),
    )

    het_samples, hom_samples = get_non_ref_samples(variant=var, samples=samples)

    # set the class attributes
    boolean_categories = [key for key in info.keys() if key.startswith('categoryboolean')]

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

    if 'blacklist' in get_config()['filter']:
        blacklist = read_json_from_path(get_config()['filter']['blacklist'])
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

    get_logger().info(f'Contig {contig} contained {contig_variants} variants')
    get_logger().info(f'Contig {contig} contained {len(contig_dict)} genes')

    return contig_dict


def read_json_from_path(
    bucket_path: str | CPGPathType | None,
    default: Any = None,
    return_model: HistoricVariants | HistoricPanels | ResultData | PanelApp | PhenotypeMatchedPanels | None = None,
) -> list | None | HistoricVariants | HistoricPanels | ResultData | PanelApp | PhenotypeMatchedPanels:
    """
    take a path to a JSON file, read into an object
    if the path doesn't exist - return the default object

    Args:
        bucket_path (str):
        default (Any):
        return_model (pydantic Models): any Model to read/validate as

    Returns:
        either the object from the JSON file, or None
    """

    if bucket_path is None:
        return default

    if isinstance(bucket_path, str):
        bucket_path = to_path(bucket_path)

    if isinstance(bucket_path, CPGPathType) and bucket_path.exists():
        with bucket_path.open() as handle:
            json_data = json.load(handle)
            if return_model:
                return return_model.model_validate(json_data)
            return json_data

    if default is not None:
        return default

    raise ValueError(f'No data found at {bucket_path}')


def get_simple_moi(input_moi: str | None, chrom: str) -> str:
    """
    takes the vast range of PanelApp MOIs, and reduces to a
    range of cases which can be easily implemented in RD analysis

    Args:
        input_moi ():
        chrom ():
    """
    if input_moi in IRRELEVANT_MOI:
        raise ValueError("unknown and other shouldn't reach this method")

    default = 'Hemi_Bi_In_Female' if chrom in X_CHROMOSOME else 'Biallelic'

    if input_moi is None:
        input_moi = default

    input_list = input_moi.translate(str.maketrans('', '', punctuation)).split()

    match input_list:
        case ['biallelic', *_additional]:
            return 'Biallelic'
        case ['both', *_additional]:
            return 'Mono_And_Biallelic'
        case ['monoallelic', *_additional]:
            if chrom in X_CHROMOSOME:
                return 'Hemi_Mono_In_Female'
            return 'Monoallelic'
        case [
            'xlinked',
            *additional,
        ] if 'biallelic' in additional:
            return 'Hemi_Bi_In_Female'
        case ['xlinked', *_additional]:
            return 'Hemi_Mono_In_Female'
        case _:
            return default


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
        if genotype_int == 1:
            het_samples.add(sam)
        if genotype_int == 3:
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
    csq_categories = get_config()['csq']['csq_string']

    # iterate over all consequences, and make each into a dict
    txc_dict = [dict(zip(csq_categories, each_csq.split('|'))) for each_csq in csq_contents.split(',')]

    # update this String to be either a float, or missing
    for each_dict in txc_dict:
        am_path = each_dict.get('am_pathogenicity')
        each_dict['am_pathogenicity'] = float(am_path) if am_path else ''

    return txc_dict


def find_comp_hets(var_list: list[VARIANT_MODELS], pedigree: peddy.Ped) -> CompHetDict:
    """
    manual implementation to find compound hets
    variants provided in the format

    [var1, var2, ..]

    generate pair content in the form
    {
        sample: {
            var_as_string: [partner_variant, ...]
        }
    }

    Args:
        var_list (list[VARIANT_MODELS]): all variants in this gene
        pedigree (): Peddy.Ped
    """

    # create an empty dictionary
    comp_het_results: CompHetDict = defaultdict(dict)  # type: ignore

    # use combinations_with_replacement to find all gene pairs
    for var_1, var_2 in combinations_with_replacement(var_list, 2):
        assert var_1.coordinates.chrom == var_2.coordinates.chrom

        if (var_1.coordinates == var_2.coordinates) or var_1.coordinates.chrom in NON_HOM_CHROM:
            continue

        # iterate over any samples with a het overlap
        for sample in var_1.het_samples.intersection(var_2.het_samples):
            phased = False
            ped_sample = pedigree.get(sample)

            # don't assess male compound hets on sex chromosomes
            if ped_sample.sex == 'male' and var_1.coordinates.chrom in X_CHROMOSOME:
                continue

            # check for both variants being in the same phase set
            if sample in var_1.phased and sample in var_2.phased:
                # check for presence of the same phase set
                for phase_set in [ps for ps in var_1.phased[sample].keys() if ps in var_2.phased[sample].keys()]:
                    if var_1.phased[sample][phase_set] == var_2.phased[sample][phase_set]:
                        phased = True
            if not phased:
                comp_het_results[sample].setdefault(var_1.coordinates.string_format, []).append(var_2)
                comp_het_results[sample].setdefault(var_2.coordinates.string_format, []).append(var_1)

    return comp_het_results


def generate_fresh_latest_results(current_results: ResultData, dataset: str, prefix: str = ''):
    """
    This will be called if a cohort has no latest results, but has
    indicated a place to save them.
    Args:
        current_results (ResultData): results from this current run
        dataset (str): dataset name for sourcing config section
        prefix (str): optional prefix for the filename
    """

    new_history = HistoricVariants(metadata=CategoryMeta(categories=get_config().get('categories', {})))
    for sample, content in current_results.results.items():
        for var in content.variants:
            new_history.results.setdefault(sample, {})[var.var_data.coordinates.string_format] = HistoricSampleVariant(
                categories={cat: var.first_seen for cat in var.categories},
                support_vars=list(var.support_vars),
                independent=var.independent,
            )
    save_new_historic(results=new_history, prefix=prefix, dataset=dataset)


def filter_results(results: ResultData, singletons: bool, dataset: str):
    """
    loads the most recent prior result set (if it exists)
    annotates previously seen variants with the most recent date seen
    write two files (total, and latest - previous)

    Args:
        results (ResultData): the results produced during this run
        singletons (bool): whether to read/write a singleton specific file
        dataset (str): dataset name for sourcing config section

    Returns: same results annotated with date-first-seen
    """

    historic_folder = get_cohort_seq_type_conf(dataset).get('historic_results')

    if historic_folder is None:
        get_logger().info('No historic data folder, no filtering')
        return

    get_logger().info('Attempting to filter current results against historic')

    # get the latest result file from the folder
    # this will be none if the folder doesn't exist or is empty
    prefix = 'singletons_' if singletons else ''

    # 2 is the required prefix, i.e. 2022_*, to discriminate vs. 'singletons_'
    # in 1000 years this might cause a problem :/ \s
    latest_results_path = find_latest_file(dataset=dataset, start=prefix or '2')

    get_logger().info(f'latest results: {latest_results_path}')

    # get latest results as a HistoricVariants object, or fail - on fail, return
    if (
        latest_results := read_json_from_path(latest_results_path, return_model=HistoricVariants)  # type: ignore
    ) is None:
        # generate and write some new latest data
        generate_fresh_latest_results(current_results=results, prefix=prefix, dataset=dataset)
        return

    assert isinstance(latest_results, HistoricVariants)

    date_annotate_results(results, latest_results)
    save_new_historic(results=latest_results, prefix=prefix, dataset=dataset)


def save_new_historic(results: HistoricVariants | HistoricPanels, dataset: str, prefix: str = ''):
    """
    save the new results in the historic results dir

    Args:
        results (): object to save as JSON
        dataset (str): the dataset to save results for
        prefix (str): name prefix for this file (optional)
    """

    directory = get_cohort_seq_type_conf(dataset).get('historic_results')
    if directory is None:
        get_logger().info('No historic data folder, no saving')
        return

    new_file = to_path(directory) / f'{prefix}{TODAY}.json'
    with new_file.open('w') as handle:
        handle.write(results.model_dump_json(indent=4))

    get_logger().info(f'Wrote new data to {new_file}')


def find_latest_file(dataset: str, results_folder: str | None = None, start: str = '', ext: str = 'json') -> str | None:
    """
    takes a directory of files, and finds the latest
    Args:
        dataset (): the dataset to fetch results for
        results_folder (): local or remote folder
        start (str): the start of the filename, if applicable
        ext (): the type of files we're looking for

    Returns: most recent file path, or None
    """

    if results_folder is None:
        results_folder = get_cohort_seq_type_conf(dataset).get('historic_results')
        if results_folder is None:
            get_logger().info('`historic_results` not present in config')
            return None

    get_logger().info(f'Using results from {results_folder}')

    date_sorted_files = sorted(
        to_path(results_folder).glob(f'{start}*.{ext}'),
        key=lambda x: x.stat().st_mtime,
        reverse=True,
    )
    if not date_sorted_files:
        return None

    return str(date_sorted_files[0].absolute())


def date_annotate_results(current: ResultData, historic: HistoricVariants):
    """
    takes the current data, and annotates with previous dates if found
    build/update the historic data within the same loop
    much simpler logic overall

    logically we can't reach this method in production without historic
    data being present. This method no longer permits historic data to be
    missing, and edits objects in place

    Args:
        current (ResultData): results generated during this run
        historic (HistoricVariants): optionally, historic data

    Returns:
        updated/instantiated cumulative data
    """

    # update to latest category descriptions
    historic.metadata.categories.update(get_config()['categories'])

    for sample, content in current.results.items():
        sample_historic = historic.results.get(sample, {})

        # check each variant found in this round
        for var in content.variants:
            var_id = var.var_data.coordinates.string_format
            current_cats = var.categories

            # this variant was previously seen
            if hist := sample_historic.get(var_id):
                # bool if this was ever independent
                # if newly independent, bump the date for current assignments
                if var.independent and not hist.independent:
                    hist.independent = True
                    for cat in current_cats:
                        hist.categories[cat] = get_granular_date()

                # if we have any new categories don't alter the date
                if new_cats := current_cats - set(hist.categories):
                    # add any new categories
                    for cat in new_cats:
                        hist.categories[cat] = get_granular_date()

                # same categories, new support
                elif new_sups := [sup for sup in var.support_vars if sup not in hist.support_vars]:
                    hist.support_vars.extend(new_sups)

                # otherwise alter the first_seen date
                else:
                    # choosing to take the latest _new_ category date
                    var.first_seen = sorted(hist.categories.values(), reverse=True)[0]

            # totally new variant
            else:
                historic.results.setdefault(sample, {})[var_id] = HistoricSampleVariant(
                    categories={cat: get_granular_date() for cat in current_cats},
                    support_vars=list(var.support_vars),
                    independent=var.independent,
                )
