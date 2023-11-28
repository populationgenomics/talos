"""
classes and methods shared across reanalysis components
"""

import time
from collections import defaultdict
from datetime import datetime
from itertools import chain, combinations_with_replacement, islice
from pathlib import Path
from string import punctuation
from typing import Any

import json
import re

import cyvcf2
import peddy
import requests

from cpg_utils import to_path, Path as CPGPathType
from cpg_utils.config import get_config

from reanalysis.models import (
    Coordinates,
    PanelApp,
    ReportVariant,
    SmallVariant,
    StructuralVariant,
    FileTypes,
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

# CONFIG_FIELDS = ['workflow']  # , 'filter', 'panels', 'categories']
# assert all(field in get_config(False).keys() for field in CONFIG_FIELDS)


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


#
# @dataclass
# class Coordinates:
#     """
#     a home for the positional variant attributes
#     """
#
#     chrom: str
#     pos: int
#     ref: str
#     alt: str
#
#     @property
#     def string_format(self) -> str:
#         """
#         forms a string representation: chr-pos-ref-alt
#         """
#         return f'{self.chrom}-{self.pos}-{self.ref}-{self.alt}'
#
#     def __lt__(self, other) -> bool:
#         """
#         enables positional sorting
#         """
#         # this will return False for same chrom and position
#         if self.chrom == other.chrom:
#             return self.pos < other.pos
#         # otherwise take the relative index from sorted chromosomes list
#         if self.chrom in CHROM_ORDER and other.chrom in CHROM_ORDER:
#             return CHROM_ORDER.index(self.chrom) < CHROM_ORDER.index(other.chrom)
#         # if self is on a canonical chromosome, sort before HLA/Decoy etc.
#         if self.chrom in CHROM_ORDER:
#             return True
#         return False
#
#     def __eq__(self, other) -> bool:
#         """
#         equivalence check
#         Args:
#             other (Coordinates):
#
#         Returns:
#             true if self == other
#
#         """
#         return (
#             self.chrom == other.chrom
#             and self.pos == other.pos
#             and self.ref == other.ref
#             and self.alt == other.alt
#         )


def get_json_response(url, max_retries=4, base_delay=1, max_delay=32):
    """
    takes a request URL, checks for healthy response, returns the JSON
    For this purpose we only expect a dictionary return
    List use-case (activities endpoint) no longer supported

    Args:
        url (str): URL to retrieve JSON format data from
        max_retries (int): maximum number of retries
        base_delay (int): initial delay between retries
        max_delay (int): maximum delay between retries

    Returns:
        the JSON response from the endpoint
    """

    retries = 0
    while retries < max_retries:
        try:
            response = requests.get(
                url, headers={'Accept': 'application/json'}, timeout=60
            )
            response.raise_for_status()  # Raise an exception for bad responses (4xx and 5xx)
            return response.json()
        except (requests.RequestException, TimeoutError) as e:
            get_logger().error(f'Request failed: {e}')
            retries += 1
            if retries < max_retries:
                delay = min(base_delay * 2**retries, max_delay)
                get_logger().warning(f'Retrying in {delay} seconds...')
                time.sleep(delay)

    raise TimeoutError('Max retries reached. Request failed.')


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


def get_new_gene_map(
    panelapp_data: PanelApp,
    pheno_panels: dict = None,
    dataset: str | None = None,
) -> dict[str, str]:
    """
    The aim here is to generate a list of all the samples for whom
    a given gene should be treated as new during this analysis. This
    prevents the need for back-filtering results at the end of
    classification.

    Generate a map of
    {gene: [samples, where, this, is, 'new']}
    """

    # find the dataset-specific panel data, if present
    # add the 'core' panel to it
    config_cohort_panels: list[int] = get_cohort_config(dataset).get(
        'cohort_panels', []
    )
    cohort_panels = config_cohort_panels + [get_config()['panels']['default_panel']]

    # collect all genes new in at least one panel
    new_genes = {
        ensg: content['new']
        for ensg, content in panelapp_data['genes'].items()
        if content['new']
    }

    # if there's no panel matching, new applies to everyone
    if pheno_panels is None:
        return {ensg: 'all' for ensg in new_genes.keys()}

    # if we have pheno-matched participants, more complex
    panel_samples = defaultdict(set)

    # double layered iteration, but only on a small object
    for sample, data in pheno_panels.items():
        for panel in data['panels']:
            panel_samples[panel].add(sample)

    pheno_matched_new = {}

    # iterate over the new genes and find out who they are new for
    for gene, panels in new_genes.items():
        if any(panel in cohort_panels for panel in panels):
            pheno_matched_new[gene] = 'all'
            continue

        # else, find the specific samples
        samples = set()
        for panel_id in panels:
            if panel_id not in panel_samples:
                raise AssertionError(f'PanelID {panel_id} not attached to any samples')
            samples.update(panel_samples[panel_id])
        pheno_matched_new[gene] = ','.join(sorted(samples))

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
        for sample, phase, genotype in zip(
            samples, map(int, var.format('PS')), var.genotypes
        ):
            # cyvcf2.Variant holds two ints, and a bool
            allele_1, allele_2, phased = genotype
            if not phased:
                continue
            gt = f'{allele_1}|{allele_2}'
            # phase set is a number
            if phase != PHASE_SET_DEFAULT:
                phased_dict[sample][phase] = gt
    except KeyError:
        get_logger().info('failed to find PS phase attributes')
        try:
            # retry using PGT & PID
            for sample, phase_gt, phase_id in zip(
                samples, var.format('PGT'), var.format('PID')
            ):
                if phase_gt != '.' and phase_id != '.':
                    phased_dict[sample][phase_id] = phase_gt
        except KeyError:
            get_logger().info('also failed using PID and PGT')

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
    new_genes: dict[str, str] | None = None,
):
    """
    takes a small variant and creates a Model from it

    Args:
        var ():
        samples ():
        as_singletons ():
        new_genes ():
    """
    coordinates = Coordinates(
        chrom=var.CHROM.replace('chr', ''), pos=var.POS, ref=var.REF, alt=var.ALT[0]
    )
    depths = dict(zip(samples, map(float, var.gt_depths)))  # type: ignore
    info: dict[str, Any] = {x.lower(): y for x, y in var.INFO} | {
        'seqr_link': coordinates.string_format
    }
    het_samples, hom_samples = get_non_ref_samples(variant=var, samples=samples)

    # hot-swap cat 2 from a boolean to a sample list - if appropriate
    if info.get('categoryboolean2', 0):
        new_gene_samples = new_genes.get(info.get('gene_id'), '')

        # if 'all', keep cohort-wide boolean flag
        if new_gene_samples == 'all':
            get_logger().debug('New applies to all samples')

        # otherwise assign only a specific sample list
        elif new_gene_samples:
            _boolcat = info.pop('categoryboolean2')
            info['categorysample2'] = new_gene_samples

        # else just remove it - shouldn't happen in prod
        else:
            _boolcat = info.pop('categoryboolean2')

    # set the class attributes
    boolean_categories = [
        key for key in info.keys() if key.startswith('categoryboolean')
    ]
    sample_categories = [key for key in info.keys() if key.startswith('categorysample')]
    sample_support = [key for key in info.keys() if key.startswith('categorysupport')]

    # overwrite with true booleans
    for cat in sample_support + boolean_categories:
        info[cat] = info.get(cat, 0) == 1

    # sample categories are a list of strings or 'missing'
    # if cohort runs as singletons, remove possibility of de novo
    # if not singletons, split each into a list of sample IDs
    for sam_cat in sample_categories:
        if as_singletons and sam_cat in REMOVE_IN_SINGLETONS:
            info[sam_cat] = []
        else:
            info[sam_cat] = (
                info[sam_cat].split(',') if info[sam_cat] != 'missing' else []
            )

    # organise PM5
    info = organise_pm5(info)
    phased = get_phase_data(samples, var)
    ab_ratios = dict(zip(samples, map(float, var.gt_alt_freqs)))
    transcript_consequences = extract_csq(csq_contents=info.pop('csq', []))

    return SmallVariant(
        coordinates=coordinates,
        info=info,
        het_samples=het_samples,
        hom_samples=hom_samples,
        boolean_categories=boolean_categories,
        sample_categories=sample_categories,
        sample_support=sample_support,
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

    # this is the right ID for Seqr
    info['seqr_link'] = info['variantid']

    coordinates = Coordinates(
        chrom=var.CHROM.replace('chr', ''),
        pos=var.POS,
        ref=var.ALT[0],
        alt=info['svlen'],
    )

    het_samples, hom_samples = get_non_ref_samples(variant=var, samples=samples)

    # set the class attributes
    boolean_categories = [
        key for key in info.keys() if key.startswith('categoryboolean')
    ]

    # overwrite with true booleans
    for cat in boolean_categories:
        info[cat] = info.get(cat, 0) == 1

    phased = get_phase_data(samples, var)

    return StructuralVariant(
        coordinates=coordinates,
        info=info,
        het_samples=het_samples,
        hom_samples=hom_samples,
        boolean_categories=boolean_categories,
        phased=phased,
    )


# CompHetDict structure: {sample: {variant_string: [variant, ...]}}
# sample: string, e,g, CGP12345
CompHetDict = dict[str, dict[str, list[SmallVariant | StructuralVariant]]]
GeneDict = dict[str, list[SmallVariant | StructuralVariant]]


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
    new_gene_map: dict[str, str],
    singletons: bool = False,
    second_source=None,
) -> GeneDict:
    """
    takes a cyvcf2.VCFReader instance, and a specified chromosome
    iterates over all variants in the region, and builds a lookup

    optionally takes a second VCF and incorporates into same dict

    Args:
        contig (): contig name from VCF header
        variant_source (): the VCF reader instance
        new_gene_map ():
        singletons ():
        second_source (): an optional second VCF (SV)

    Returns:
        A lookup in the form
        {
            gene1: [var1, var2],
            gene2: [var3],
            ...
        }
    """

    if 'blacklist' in get_config()['filter']:
        blacklist = read_json_from_path(get_config()['filter']['blacklist'])
    else:
        blacklist = []

    # a dict to allow lookup of variants on this whole chromosome
    contig_variants = 0
    contig_dict = defaultdict(list)

    # iterate over all variants on this contig and store by unique key
    # if contig has no variants, prints an error and returns []
    for variant in variant_source(contig):

        abs_var = create_small_variant(
            var=variant,
            samples=variant_source.samples,
            as_singletons=singletons,
            new_genes=new_gene_map,
        )

        if abs_var.coordinates.string_format in blacklist:
            get_logger().info(
                f'Skipping blacklisted variant: {abs_var.coordinates.string_format}'
            )
            continue

        # if unclassified, skip the whole variant
        if not abs_var.is_classified:
            continue

        # update the variant count
        contig_variants += 1

        # update the gene index dictionary
        contig_dict[abs_var.info.get('gene_id')].append(abs_var)

    if second_source:
        second_source_variants = 0
        for variant in second_source(contig):
            # create an abstract SV variant
            abs_var = create_structural_variant(
                var=variant, samples=second_source.samples
            )
            # update the variant count
            second_source_variants += 1

            # update the gene index dictionary
            contig_dict[abs_var.info.get('gene_id')].append(abs_var)

        get_logger().info(
            f'Contig {contig} contained {second_source_variants} variants'
        )

    get_logger().info(f'Contig {contig} contained {contig_variants} variants')
    get_logger().info(f'Contig {contig} contained {len(contig_dict)} genes')

    return contig_dict


def read_json_from_path(
    bucket_path: str | CPGPathType | None, default: Any = None
) -> dict | list | None:
    """
    take a path to a JSON file, read into an object
    if the path doesn't exist - return the default object

    Args:
        bucket_path (str):
        default (Any):

    Returns:
        either the object from the JSON file, or None
    """

    if bucket_path is None:
        return default

    if isinstance(bucket_path, str):
        bucket_path = to_path(bucket_path)

    if isinstance(bucket_path, CPGPathType) and bucket_path.exists():
        with bucket_path.open() as handle:
            return json.load(handle)
    return default


def write_output_json(output_path: str, object_to_write: dict):
    """
    writes object to a json file, to_path provides platform abstraction

    Args:
        output_path ():
        object_to_write ():
    """

    get_logger().info(f'Writing output JSON file to {output_path}')
    out_route = to_path(output_path)

    if out_route.exists():
        get_logger().info(f'Output path {output_path!r} exists, will be overwritten')

    with out_route.open('w') as fh:
        json.dump(object_to_write, fh, indent=4, default=list)


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


def extract_csq(csq_contents) -> list[dict]:
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
    txc_dict = [
        dict(zip(csq_categories, each_csq.split('|')))
        for each_csq in csq_contents.split(',')
    ]

    # update this String to be either a float, or missing
    for each_dict in txc_dict:
        am_path = each_dict.get('am_pathogenicity')
        each_dict['am_pathogenicity'] = float(am_path) if am_path else ''

    return txc_dict


def find_comp_hets(
    var_list: list[SmallVariant | StructuralVariant], pedigree: peddy.Ped
) -> CompHetDict:
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
        var_list (list[SmallVariant | StructuralVariant]): all variants in this gene
        pedigree (): Peddy.Ped
    """

    # create an empty dictionary
    comp_het_results: CompHetDict = defaultdict(dict)  # type: ignore

    # use combinations_with_replacement to find all gene pairs
    for var_1, var_2 in combinations_with_replacement(var_list, 2):
        assert var_1.coordinates.chrom == var_2.coordinates.chrom

        if (
            var_1.coordinates == var_2.coordinates
        ) or var_1.coordinates.chrom in NON_HOM_CHROM:
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
                for phase_set in [
                    ps
                    for ps in var_1.phased[sample].keys()
                    if ps in var_2.phased[sample].keys()
                ]:
                    if (
                        var_1.phased[sample][phase_set]
                        == var_2.phased[sample][phase_set]
                    ):
                        phased = True
            if not phased:
                comp_het_results[sample].setdefault(
                    var_1.coordinates.string_format, []
                ).append(var_2)
                comp_het_results[sample].setdefault(
                    var_2.coordinates.string_format, []
                ).append(var_1)

    return comp_het_results


def filter_results(results: dict, singletons: bool, dataset: str) -> dict:
    """
    loads the most recent prior result set (if it exists)
    annotates previously seen variants with the most recent date seen
    write two files (total, and latest - previous)

    Args:
        results (): the results produced during this run
        singletons (bool): whether to read/write a singleton specific file
        dataset (str): dataset name for sourcing config section

    Returns: same results annotated with date-first-seen
    """

    historic_folder = get_cohort_seq_type_conf(dataset).get('historic_results')

    if historic_folder is None:
        get_logger().info('No historic data folder, no filtering')
        # results, _cumulative = date_annotate_results(results)
        return results

    get_logger().info('Attempting to filter current results against historic')

    # get the latest result file from the folder
    # this will be none if the folder doesn't exist or is empty
    prefix = 'singletons_' if singletons else ''

    # 2 is the required prefix, i.e. 2022_*, to discriminate vs. 'singletons_'
    # in 1000 years this might cause a problem :/ \s
    latest_results_path = find_latest_file(dataset=dataset, start=prefix or '2')

    get_logger().info(f'latest results: {latest_results_path}')

    latest_results: dict = read_json_from_path(latest_results_path)  # type: ignore

    results, cumulative = date_annotate_results(results, latest_results)
    save_new_historic(results=cumulative, prefix=prefix, dataset=dataset)

    return results


def save_new_historic(
    results: dict, dataset: str, prefix: str = '', directory: str | None = None
):
    """
    save the new results in the historic results dir

    Args:
        results (): object to save as a JSON file
        dataset (str): the dataset to save results for
        prefix (str): name prefix for this file (optional)
        directory (): defaults to historic_data from config
    """

    if directory is None:
        directory = get_cohort_seq_type_conf(dataset).get('historic_results')
        if directory is None:
            get_logger().info('No historic results directory, nothing written')
            return

    new_file = to_path(directory) / f'{prefix}{TODAY}.json'
    with new_file.open('w') as handle:
        json.dump(results, handle, indent=4, default=list)

    get_logger().info(f'Wrote new data to {new_file}')


def find_latest_file(
    dataset: str, results_folder: str | None = None, start: str = '', ext: str = 'json'
) -> str | None:
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


def date_annotate_results(
    current: dict[str, dict | list[ReportVariant]], historic: dict | None = None
) -> tuple[dict, dict]:
    """
    takes the current data, and annotates with previous dates if found
    build/update the historic data within the same loop
    much simpler logic overall

    Args:
        current (dict): results generated during this run
        historic (): optionally, historic data

    Returns: date-annotated results and cumulative data
    """

    # if there's no historic data, make some
    if historic is None:
        historic = {
            'metadata': {'categories': get_config()['categories']},
            'results': {},
        }

    # update to latest format
    elif 'results' not in historic.keys():
        historic = {
            'metadata': {'categories': get_config()['categories']},
            'results': historic,
        }

    # update to latest category descriptions
    historic['metadata'].setdefault('categories', {}).update(get_config()['categories'])

    for sample, content in current.items():

        # totally absent? start populating for this sample
        if sample not in historic['results']:
            historic['results'][sample] = {}

        # check each variant found in this round
        for var in content['variants']:  # type: ignore
            var_id = var.var_data.coordinates.string_format
            current_cats = set(var.var_data.categories)

            # this variant was previously seen
            if var_id in historic['results'][sample]:

                hist = historic['results'][sample][var_id]

                # bool if this was ever independent
                hist['independent'] = var.independent or hist.get('independent', False)

                # if we have any new categories don't alter the date
                if new_cats := current_cats - set(hist['categories'].keys()):

                    # add any new categories
                    for cat in new_cats:
                        hist['categories'][cat] = get_granular_date()

                # same categories, new support
                elif new_sups := [
                    sup for sup in var.support_vars if sup not in hist['support_vars']
                ]:
                    hist['support_vars'].extend(new_sups)

                # otherwise alter the first_seen date
                else:
                    # choosing to take the latest _new_ category date
                    recent = sorted(hist['categories'].values(), reverse=True)[0]
                    var.first_seen = recent

            # totally new variant
            else:
                historic['results'][sample][var_id] = {
                    'categories': {cat: get_granular_date() for cat in current_cats},
                    'support_vars': var.support_vars,
                    'independent': var.independent,
                }

    return current, historic
