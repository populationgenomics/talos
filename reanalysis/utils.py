"""
classes and methods shared across reanalysis components
"""


from collections import defaultdict
from dataclasses import dataclass, is_dataclass, field
from datetime import datetime
from enum import Enum
from itertools import chain, combinations_with_replacement, islice
from pathlib import Path
from string import punctuation
from typing import Any

import json
import logging
import re
import os

import requests

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.git import get_git_repo_root

# pylint: disable=too-many-lines


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
GRANULAR = datetime.now().strftime('%Y-%m-%d')

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


class FileTypes(Enum):
    """
    enumeration of permitted input file types
    """

    HAIL_TABLE = '.ht'
    MATRIX_TABLE = '.mt'
    VCF = '.vcf'
    VCF_GZ = '.vcf.gz'
    VCF_BGZ = '.vcf.bgz'


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

@dataclass
class Coordinates:
    """
    a home for the positional variant attributes
    """

    chrom: str
    pos: int
    ref: str
    alt: str

    @property
    def string_format(self):
        """
        forms a string representation
        chr-pos-ref-alt
        """
        return f'{self.chrom}-{self.pos}-{self.ref}-{self.alt}'

    def __lt__(self, other):
        """
        enables positional sorting
        """
        # this will return False for same chrom and position
        if self.chrom == other.chrom:
            return self.pos < other.pos
        # otherwise take the relative index from sorted chromosomes list
        if self.chrom in CHROM_ORDER and other.chrom in CHROM_ORDER:
            return CHROM_ORDER.index(self.chrom) < CHROM_ORDER.index(other.chrom)
        # if self is on a canonical chromosome, sort before HLA/Decoy etc.
        if self.chrom in CHROM_ORDER:
            return True
        return False

    def __eq__(self, other):
        """
        equivalence check
        Args:
            other (Coordinates):

        Returns:
            true if self == other

        """
        return (
            self.chrom == other.chrom
            and self.pos == other.pos
            and self.ref == other.ref
            and self.alt == other.alt
        )


def get_json_response(url: str) -> Any:
    """
    takes a request URL, checks for healthy response, returns the JSON
    For this purpose we only expect a dictionary return
    List use-case (activities endpoint) no longer supported

    Args:
        url ():

    Returns:
        the JSON response from the endpoint
    """

    response = requests.get(url, headers={'Accept': 'application/json'}, timeout=60)
    response.raise_for_status()
    return response.json()


def get_new_gene_map(
    panelapp_data: dict, pheno_panels: dict | None = None
) -> dict[str, str]:
    """
    The aim here is to generate a list of all the samples for whom
    a given gene should be treated as new during this analysis. This
    prevents the need for back-filtering results at the end of
    classification.

    Generate a map of
    { gene: [samples, where, this, is, 'new']}
    """

    # pull out the core panel once
    core_panel = get_config()['workflow']['default_panel']

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
        if core_panel in panels:
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
    phased_dict = defaultdict(dict)

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
        logging.info('failed to find PS phase attributes')
        try:
            # retry using PGT & PID
            for sample, phase_gt, phase_id in zip(
                samples, var.format('PGT'), var.format('PID')
            ):
                if phase_gt != '.' and phase_id != '.':
                    phased_dict[sample][phase_id] = phase_gt
        except KeyError:
            logging.info('also failed using PID and PGT')

    return dict(phased_dict)


@dataclass
class AbstractVariant:  # pylint: disable=too-many-instance-attributes
    """
    create class to contain all content from cyvcf2 object
    """

    def __init__(
        self,
        var,
        samples: list[str],
        as_singletons=False,
        new_genes: dict[str, str] | None = None,
    ):
        """
        Args:
            var (cyvcf2.Variant):
            samples (list):
            as_singletons (bool):
            new_genes (dict):
        """

        # extract the coordinates into a separate object
        self.coords = Coordinates(
            var.CHROM.replace('chr', ''), var.POS, var.REF, var.ALT[0]
        )

        # get all zygosities once per variant
        # abstraction avoids pulling per-sample calls again later
        self.het_samples, self.hom_samples = get_non_ref_samples(
            variant=var, samples=samples
        )

        # overwrite the non-standard cyvcf2 representation
        self.info: dict[str, Any] = {x.lower(): y for x, y in var.INFO}

        # hot-swap cat 2 from a boolean to a sample list - if appropriate
        if self.info.get('categoryboolean2', 0):
            new_gene_samples = new_genes.get(self.info.get('gene_id'), '')

            # if 'all', keep cohort-wide boolean flag
            if new_gene_samples == 'all':
                logging.debug('New applies to all samples')

            # otherwise assign only a specific sample list
            elif new_gene_samples:
                _boolcat = self.info.pop('categoryboolean2')
                self.info['categorysample2'] = new_gene_samples

            # else just remove it - shouldn't happen in prod
            else:
                _boolcat = self.info.pop('categoryboolean2')

        # set the class attributes
        self.boolean_categories = [
            key for key in self.info.keys() if key.startswith('categoryboolean')
        ]
        self.sample_categories = [
            key for key in self.info.keys() if key.startswith('categorysample')
        ]
        self.sample_support = [
            key for key in self.info.keys() if key.startswith('categorysupport')
        ]

        # overwrite with true booleans
        for cat in self.sample_support + self.boolean_categories:
            self.info[cat] = self.info.get(cat, 0) == 1

        # sample categories are a list of strings or 'missing'
        # if cohort runs as singletons, remove possibility of de novo
        # if not singletons, split each into a list of sample IDs
        for sam_cat in self.sample_categories:
            if as_singletons and sam_cat in REMOVE_IN_SINGLETONS:
                self.info[sam_cat] = []
            else:
                self.info[sam_cat] = (
                    self.info[sam_cat].split(',')
                    if self.info[sam_cat] != 'missing'
                    else []
                )

        self.transcript_consequences: list[dict[str, str]] = []
        if 'csq' in self.info:
            self.transcript_consequences = extract_csq(
                csq_contents=self.info.pop('csq')
            )

        # identify variant sets phased with this one
        # cyvcf2 uses a default value for the phase set, skip that
        # this is restricted to a single int for phase_set
        try:
            self.phased = get_phase_data(samples, var)
        except KeyError:
            self.phased = {} 

        self.ab_ratios = dict(zip(samples, map(float, var.gt_alt_freqs)))
        self.categories = []

    def __str__(self):
        return repr(self)

    def __lt__(self, other):
        return self.coords < other.coords

    def __eq__(self, other):
        return self.coords == other.coords

    @property
    def has_boolean_categories(self) -> bool:
        """
        check that the variant has at least one assigned class
        """
        return any(self.info[value] for value in self.boolean_categories)

    @property
    def has_sample_categories(self) -> bool:
        """
        check that the variant has any list-category entries
        """
        return any(self.info[value] for value in self.sample_categories)

    @property
    def has_support(self) -> bool:
        """
        check for a True flag in any CategorySupport* attribute
        Returns:
            True if variant is support
        """
        return any(self.info[value] for value in self.sample_support)

    @property
    def category_non_support(self) -> bool:
        """
        check the variant has at least one non-support category assigned
        Returns:
            True if var has a non-support category assigned
        """
        return self.has_sample_categories or self.has_boolean_categories

    @property
    def is_classified(self) -> bool:
        """
        check for at least one assigned class, inc. support
        Returns:
            True if classified
        """
        return self.category_non_support or self.has_support

    @property
    def support_only(self) -> bool:
        """
        check that the variant is exclusively cat. support
        Returns:
            True if support only
        """
        return self.has_support and not self.category_non_support

    def category_values(self, sample: str) -> list[str]:
        """
        get all variant categories; sample-specific checks for de novo

        Args:
            sample (str): sample id

        Returns:
            list of all categories applied to this variant
        """
        categories = [
            bool_cat.replace('categoryboolean', '')
            for bool_cat in self.boolean_categories
            if self.info[bool_cat]
        ]
        if self.has_support:
            categories.append('support')

        if self.sample_de_novo(sample_id=sample):
            categories.append('de_novo')

        # mutually exlusive with the boolean category2 value
        if new := self.info.get('categorysample2'):
            if any(x in new for x in ['all', sample]):
                categories.append('2')

        return categories

    def sample_de_novo(self, sample_id: str) -> bool:
        """
        check if variant is de novo for this sample

        Args:
            sample_id (str):

        Returns:
            bool: True if this sample forms de novo
        """
        return sample_id in self.info.get('categorysample4', [])

    def sample_categorised_check(self, sample_id: str) -> bool:
        """
        check if any *sample categories applied for this sample

        Args:
            sample_id (str):

        Returns:
            bool: True if this sample features in any
                  named-sample category, includes 'all'
        """
        return any(
            sam in self.info[sam_cat]
            for sam_cat in self.sample_categories
            for sam in [sample_id, 'all']
        )

    def sample_category_check(self, sample_id: str, allow_support: bool = True) -> bool:
        """
        take a specific sample and check for assigned categories
        optionally, include checks for support category

        Args:
            sample_id (str):
            allow_support: (bool) also check for support

        Returns:
            True if the variant is categorised for this sample
        """
        big_cat = self.category_non_support or self.sample_categorised_check(sample_id)
        if allow_support:
            return big_cat or self.has_support
        return big_cat

    def get_sample_flags(self, sample: str) -> list[str]:
        """
        gets all report flags for this sample - currently only one flag
        """
        return self.check_ab_ratio(sample)

    def check_ab_ratio(self, sample: str) -> list[str]:
        """
        AB ratio test for this sample's variant call
        Args:
            sample (str): sample ID

        Returns:
            list[str]: empty, or indicating an AB ratio failure
        """
        het = sample in self.het_samples
        hom = sample in self.hom_samples
        variant_ab = self.ab_ratios.get(sample, 0.0)
        if (
            (variant_ab <= 0.15)
            or (het and not 0.25 <= variant_ab <= 0.75)
            or (hom and variant_ab <= 0.85)
        ):
            return ['AB Ratio']
        return []


# CompHetDict structure: {sample: {variant_string: [variant, ...]}}
# sample: string, e,g, CGP12345
CompHetDict = dict[str, dict[str, list[AbstractVariant]]]
GeneDict = dict[str, list[AbstractVariant]]


@dataclass
class ReportedVariant:
    """
    minimal model representing variant categorisation event
    the initial variant
    the MOI passed
    the support (if any)
    allows for the presence of flags e.g. Borderline AB ratio
    """

    # pylint: disable=too-many-instance-attributes

    sample: str
    family: str
    gene: str
    var_data: AbstractVariant
    reasons: set[str]
    genotypes: dict[str, str]
    supported: bool = field(default=False)
    support_vars: list[str] = field(default_factory=list)
    flags: list[str] = field(default_factory=list)
    phenotypes: list[str] = field(default_factory=list)
    first_seen: str = GRANULAR

    def __eq__(self, other):
        """
        makes reported variants comparable
        """
        self_supvar = set(self.support_vars)
        other_supvar = set(other.support_vars)
        return (
            self.sample == other.sample
            and self.var_data.coords == other.var_data.coords
            and self.supported == other.supported
            and self_supvar == other_supvar
        )

    def __lt__(self, other):
        return self.var_data.coords < other.var_data.coords


def canonical_contigs_from_vcf(reader) -> set[str]:
    """
    reader is a cyvcf2.VCFReader
    read the header fields from the VCF handle
    return a set of all 'canonical' contigs
    :param reader:
    :return:
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
) -> GeneDict:
    """
    takes a cyvcf2.VCFReader instance, and a specified chromosome
    iterates over all variants in the region, and builds a lookup

    Args:
        contig (): contig name from VCF header
        variant_source (): the VCF reader instance
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

    blacklist = []
    if 'blacklist' in get_config()['filter'].keys():
        blacklist = read_json_from_path(get_config()['filter']['blacklist'])

    # a dict to allow lookup of variants on this whole chromosome
    contig_variants = 0
    contig_dict = defaultdict(list)

    # iterate over all variants on this contig and store by unique key
    # if contig has no variants, prints an error and returns []
    for variant in variant_source(contig):
        abs_var = AbstractVariant(
            var=variant,
            samples=variant_source.samples,
            as_singletons=singletons,
            new_genes=new_gene_map,
        )

        if abs_var.coords.string_format in blacklist:
            logging.info(
                f'Skipping blacklisted variant: {abs_var.coords.string_format}'
            )
            continue

        # if unclassified, skip the whole variant
        if not abs_var.is_classified:
            continue

        # update the variant count
        contig_variants += 1

        # update the gene index dictionary
        contig_dict[abs_var.info.get('gene_id')].append(abs_var)

    logging.info(f'Contig {contig} contained {contig_variants} variants')
    logging.info(f'Contig {contig} contained {len(contig_dict)} genes')
    return contig_dict


def read_json_from_path(bucket_path: str | Path | None, default: Any = None) -> Any:
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

    if bucket_path.exists():
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

    logging.info(f'Writing output JSON file to {output_path}')
    out_route = to_path(output_path)

    if out_route.exists():
        logging.info(f'Output path {output_path!r} exists, will be overwritten')

    with out_route.open('w') as fh:
        json.dump(object_to_write, fh, indent=4, default=list)


def get_simple_moi(input_moi: str | None, chrom: str) -> str | None:
    """
    takes the vast range of PanelApp MOIs, and reduces to a
    range of cases which can be easily implemented in RD analysis

    Args:
        input_moi ():
        chrom ():

    Returns:

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
            return 'Monoallelic'
        case [
            'x-linked',
            *additional,
        ] if 'biallelic' in additional:  # pylint: disable='used-before-assignment'
            return 'Hemi_Bi_In_Female'
        case ['x-linked', *_additional]:
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
    handle extraction of the CSQ entries
    :param csq_contents:
    :return:
    """

    # allow for no CSQ data, i.e. splice variant
    if not csq_contents:
        return []

    # break mono-CSQ-string into components
    csq_categories = get_config()['csq']['csq_string']

    # iterate over all consequences, and make each into a dict
    return [
        dict(zip(csq_categories, each_csq.split('|')))
        for each_csq in csq_contents.split(',')
    ]


class CustomEncoder(json.JSONEncoder):
    """
    to be used as a JSON encoding class
    - replaces all sets with lists
    - replaces dataclass objects with a dictionary of the same
    """

    def default(self, o):
        """
        takes an arbitrary object, and forms a JSON representation
        where the object doesn't have an easy string representation,
        transform to a valid object: set -> list, class -> dict
        :param o:
        :return:
        """

        if is_dataclass(o):
            return o.__dict__
        if isinstance(o, set):
            return list(o)
        return json.JSONEncoder.default(self, o)


def find_comp_hets(var_list: list[AbstractVariant], pedigree) -> CompHetDict:
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

    :param var_list:
    :param pedigree: peddy.Ped
    """

    # create an empty dictionary
    comp_het_results = defaultdict(dict)

    # use combinations_with_replacement to find all gene pairs
    for var_1, var_2 in combinations_with_replacement(var_list, 2):
        assert var_1.coords.chrom == var_2.coords.chrom

        if (var_1.coords == var_2.coords) or var_1.coords.chrom in NON_HOM_CHROM:
            continue

        # iterate over any samples with a het overlap
        for sample in var_1.het_samples.intersection(var_2.het_samples):
            phased = False
            ped_sample = pedigree.get(sample)

            # don't assess male compound hets on sex chromosomes
            if ped_sample.sex == 'male' and var_1.coords.chrom in X_CHROMOSOME:
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
                    var_1.coords.string_format, []
                ).append(var_2)
                comp_het_results[sample].setdefault(
                    var_2.coords.string_format, []
                ).append(var_1)

    return comp_het_results


def filter_results(results: dict, singletons: bool) -> dict:
    """
    loads the most recent prior result set (if it exists)
    annotates previously seen variants with the most recent date seen
    write two files (total, and latest - previous)

    Args:
        results (): the results produced during this run
        singletons (bool): whether to read/write a singleton specific file

    Returns:
        the same results back-filtered to remove previous results
    """

    historic_folder = get_config()['dataset_specific'].get('historic_results')

    if historic_folder is None:
        logging.info('No historic data folder, no filtering')
        return results

    logging.info('Attempting to filter current results against historic')

    # get the latest result file from the folder
    # this will be none if the folder doesn't exist or is empty
    prefix = 'singletons_' if singletons else ''

    # 2 is the required prefix, i.e. 2022_*, to discriminate vs. 'singletons_'
    # in 1000 years this might cause a problem :/ \s
    latest_results = find_latest_file(start=prefix or '2')

    logging.info(f'latest results: {latest_results}')

    results, cumulative = date_annotate_results(
        results, read_json_from_path(latest_results)
    )
    save_new_historic(results=cumulative, prefix=prefix)

    return results


def save_new_historic(results: dict, prefix: str = '', directory: str | None = None):
    """
    save the new results in the historic results dir

    Args:
        results (): object to save as a JSON file
        prefix (str): name prefix for this file (optional)
        directory (): defaults to historic_data from config
    """

    if directory is None:
        directory = get_config()['dataset_specific'].get('historic_results')
        if directory is None:
            logging.info('No historic results directory, nothing written')
            return

    new_file = to_path(directory) / f'{prefix}{TODAY}.json'
    with new_file.open('w') as handle:
        json.dump(results, handle, indent=4, default=list)

    logging.info(f'Wrote new data to {new_file}')


def find_latest_file(
    results_folder: str | None = None, start: str = '', ext: str = 'json'
) -> str | None:
    """
    takes a directory of files, and finds the latest
    Args:
        results_folder (): local or remote folder
        start (str): the start of the filename, if applicable
        ext (): the type of files we're looking for

    Returns:
        most recent file path, or None
    """

    if results_folder is None:
        results_folder = (
            get_config().get('dataset_specific', {}).get('historic_results')
        )
        if results_folder is None:
            logging.info('`historic_results` not present in config')
            return None

    logging.info(f'Using results from {results_folder}')

    date_sorted_files = sorted(
        to_path(results_folder).glob(f'{start}*.{ext}'),
        key=lambda x: x.stat().st_mtime,
        reverse=True,
    )
    if not date_sorted_files:
        return None

    return str(date_sorted_files[0].absolute())


def date_annotate_results(
    current: dict[str, dict | list[ReportedVariant]], historic: dict | None = None
) -> tuple[dict, dict]:
    """
    takes the current data, and annotates with previous dates if found
    build/update the historic data within the same loop
    much simpler logic overall

    Args:
        current ():
        historic (): optionally, historic data

    Returns:
        the date-annotated results and cumulative data
    """

    # if there's no historic data, make some
    if historic is None:
        historic = {}

    for sample, content in current.items():

        # totally absent? start populating for this sample
        if sample not in historic:
            historic[sample] = {}

        # check each variant found in this round
        for var in content['variants']:
            var_id = var.var_data.coords.string_format
            current_cats = set(var.var_data.categories)

            # this variant was previously seen
            if var_id in historic[sample]:

                hist = historic[sample][var_id]
                historic_cats = set(hist['categories'].keys())

                # if we have any new categories don't alter the date
                if new_cats := current_cats - historic_cats:

                    # add any new categories
                    for cat in new_cats:
                        hist['categories'][cat] = GRANULAR

                # same categories, new support
                elif new_sups := [
                    sup for sup in var.support_vars if sup not in hist['support_vars']
                ]:
                    hist['support_vars'].extend(new_sups)

                # otherwise alter the first_seen date
                # todo first_seen is the wrong nomenclature here
                else:
                    # choosing to take the latest _new_ category date
                    recent = sorted(hist['categories'].values(), reverse=True)[0]
                    var.first_seen = recent

            # totally new variant
            else:
                historic[sample][var_id] = {
                    'categories': {cat: GRANULAR for cat in current_cats},
                    'support_vars': var.support_vars,
                }

    return current, historic
