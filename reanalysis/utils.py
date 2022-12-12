"""
a collection of classes and methods
which may be shared across reanalysis components
"""

from collections import defaultdict
from dataclasses import dataclass, is_dataclass, field
from datetime import datetime
from enum import Enum
from itertools import combinations_with_replacement
from pathlib import Path
from typing import Any

import json
import logging
import re
import requests

from cpg_utils import to_path
from cpg_utils.config import get_config


HOMREF: int = 0
HETALT: int = 1
UNKNOWN: int = 2
HOMALT: int = 3

# in cyVCF2, these ints represent HOMREF, and UNKNOWN
BAD_GENOTYPES: set[int] = {HOMREF, UNKNOWN}

PHASE_SET_DEFAULT = -2147483648
CHROM_ORDER = list(map(str, range(1, 23))) + [
    'X',
    'Y',
    'MT',
    'M',
]

X_CHROMOSOME = {'X'}
NON_HOM_CHROM = {'Y', 'MT', 'M'}

TODAY = datetime.now().strftime('%Y-%m-%d_%H:%M')


class FileTypes(Enum):
    """
    enumeration of permitted input file types
    """

    HAIL_TABLE = '.ht'
    MATRIX_TABLE = '.mt'
    VCF = '.vcf'
    VCF_GZ = '.vcf.gz'
    VCF_BGZ = '.vcf.bgz'


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
    raise Exception(f'File cannot be definitively typed: {str(extensions)}')


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
        positional sorting
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


def get_json_response(url: str) -> dict[str, Any]:
    """
    takes a request URL, checks for healthy response, returns the JSON
    For this purpose we only expect a dictionary return
    List use-case (activities endpoint) no longer supported

    Args:
        url ():

    Returns:

    """

    response = requests.get(url, headers={'Accept': 'application/json'}, timeout=60)
    response.raise_for_status()
    return response.json()


def get_phase_data(samples, var) -> dict[str, dict[int, str]]:
    """
    read phase data from this variant

    Args:
        samples ():
        var ():

    Returns:

    """
    phased_dict = defaultdict(dict)

    # first set the numpy.ndarray to be a list of ints
    # the zip against ordered sample IDs
    # this might need to store the exact genotype too
    # i.e. 0|1 and 1|0 can be in the same phase-set
    # but are un-phased variants
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

    return dict(phased_dict)


@dataclass
class AbstractVariant:  # pylint: disable=too-many-instance-attributes
    """
    create a bespoke variant class
    pull all content out of the cyvcf2 object
    we could have a separate implementation for pyvcf,
    or for a direct parser...
    """

    def __init__(
        self,
        var,
        samples: list[str],
        as_singletons=False,
    ):
        """

        Args:
            var (cyvcf2.Variant):
            samples (list):
            as_singletons (bool):
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

        # de novo categories are a list of strings or 'missing'
        # if cohort runs as singletons, remove possibility of de novo
        # if not singletons, split each into a list of sample IDs
        for sam_cat in self.sample_categories:
            if as_singletons:
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
        self.phased = get_phase_data(samples, var)

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
        :return:
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
        """
        return any(self.info[value] for value in self.sample_support)

    @property
    def category_non_support(self) -> bool:
        """
        check the variant has at least one non-support category assigned
        :return:
        """
        return self.has_sample_categories or self.has_boolean_categories

    @property
    def is_classified(self) -> bool:
        """
        check that the variant has at least one assigned class
        supporting category is considered here
        :return:
        """
        return self.category_non_support or self.has_support

    @property
    def support_only(self) -> bool:
        """
        checks that the variant was only class 4
        :return:
        """
        return self.has_support and not self.category_non_support

    def category_values(self, sample: str) -> list[str]:
        """
        get a list values representing the classes present on this variant
        for each category, append that number if the class is present
        - de novo on a per-sample basis
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

        return categories

    def sample_de_novo(self, sample_id: str) -> bool:
        """
        takes a specific sample ID, to check if the sample has a de novo call

        :param sample_id:
        :return:
        """
        return any(
            sample_id in self.info[sam_cat] for sam_cat in self.sample_categories
        )

    def sample_specific_category_check(self, sample_id: str) -> bool:
        """

        :param sample_id:
        :return:
        """
        return self.category_non_support or self.sample_de_novo(sample_id)

    def get_sample_flags(self, sample: str) -> list[str]:
        """
        gets all report flags for this sample
        """
        return self.check_ab_ratio(sample)

    def check_ab_ratio(self, sample: str) -> list[str]:
        """
        AB ratio test for this sample's variant call

        ab = matrix.AD[1] / hl.sum(matrix.AD)
        return matrix.filter_entries(
            (matrix.GT.is_hom_ref() & (ab <= 0.15))
            | (matrix.GT.is_het() & (ab >= 0.25) & (ab <= 0.75))
            | (matrix.GT.is_hom_var() & (ab >= 0.85))
        )
        :param sample: this affected individual
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

    sample: str
    gene: str
    var_data: AbstractVariant
    reasons: set[str]
    supported: bool = field(default=False)
    support_vars: list[str] = field(default_factory=list)
    flags: list[str] = field(default_factory=list)

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
    contig: str, variant_source, panelapp_data: dict, singletons: bool = False
) -> GeneDict:
    """
    takes a cyvcf2.VCFReader instance, and a specified chromosome
    iterates over all variants in the region, and builds a lookup
    {
        gene1: [var1, var2],
        gene2: [var3],
        ...
    }
    :param contig: contig name from header (canonical_contigs_from_vcf)
    :param variant_source: VCF reader instance
    :param panelapp_data:
    :param singletons:
    :return: populated lookup dict
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
        )

        if abs_var.coords.string_format in blacklist:
            logging.info(
                f'Skipping blacklisted variant: {abs_var.coords.string_format}'
            )
            continue

        # if the gene isn't 'new' in any panel, remove Cat2 flag
        if abs_var.info.get('categoryboolean2'):
            abs_var.info['categoryboolean2'] = (
                len(panelapp_data[abs_var.info.get('gene_id')].get('new', [])) > 0
            )

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


def read_json_from_path(bucket_path: str) -> dict[str, Any]:
    """
    take a path to a JSON file, read into an object
    :param bucket_path:
    """
    with to_path(bucket_path).open() as handle:
        return json.load(handle)


# most lenient to most conservative
# usage = if we have two MOIs for the same gene, take the broadest
ORDERED_MOIS = [
    'Mono_And_Biallelic',
    'Monoallelic',
    'Hemi_Mono_In_Female',
    'Hemi_Bi_In_Female',
    'Biallelic',
]


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
        logging.info(f'Output path "{output_path}" exists, will be overwritten')

    with out_route.open('w') as fh:
        json.dump(object_to_write, fh, indent=4, default=str)


def get_simple_moi(panel_app_moi: str | None) -> str:
    """
    takes the vast range of PanelApp MOIs, and reduces to a reduced
    range of cases which can be easily implemented in RD analysis
    This is required to reduce the complexity of an MVP
    Could become a strict enumeration

    {
        Biallelic: [all, panelapp, moi, equal, to, biallelic],
        Monoallelic: [all MOI to be interpreted as monoallelic]
    }

    :param panel_app_moi: full PanelApp string or None
    :return: a simplified representation
    """

    # default to considering both. NOTE! Many genes have Unknown MOI!
    simple_moi = 'Mono_And_Biallelic'
    # try-except permits the moi to be None
    try:
        lower_moi = panel_app_moi.lower()
        if lower_moi is None or lower_moi == 'unknown':
            # exit iteration, all simple moi considered
            return simple_moi
    except AttributeError:
        return simple_moi

    # ideal for match-case, coming to a python 3.10 near you!
    if lower_moi.startswith('biallelic'):
        simple_moi = 'Biallelic'
    elif lower_moi.startswith('both'):
        simple_moi = 'Mono_And_Biallelic'
    elif lower_moi.startswith('mono'):
        simple_moi = 'Monoallelic'
    elif lower_moi.startswith('x-linked'):
        if 'biallelic' in panel_app_moi:
            simple_moi = 'Hemi_Bi_In_Female'
        else:
            simple_moi = 'Hemi_Mono_In_Female'

    return simple_moi


def get_non_ref_samples(variant, samples: list[str]) -> tuple[set[str], set[str]]:
    """
    variant is a cyvcf2.Variant
    for this variant, find all samples with a call
    cyvcf2 uses 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
    return het, hom, and the union of het and hom
    :param variant:
    :param samples:
    :return:
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
    :return:
    """

    # create an empty dictionary
    comp_het_results = defaultdict(dict)

    # use combinations_with_replacement to find all gene pairs
    for var_1, var_2 in combinations_with_replacement(var_list, 2):
        assert var_1.coords.chrom == var_2.coords.chrom

        if (var_1.coords == var_2.coords) or var_1.coords.chrom in NON_HOM_CHROM:
            continue

        sex_chrom = var_1.coords.chrom in X_CHROMOSOME

        # iterate over any samples with a het overlap
        for sample in var_1.het_samples.intersection(var_2.het_samples):
            phased = False
            ped_sample = pedigree.get(sample)

            # don't assess male compound hets on sex chromosomes
            if ped_sample.sex == 'male' and sex_chrom:
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
    takes a set of results
    loads the most recent prior result set (if it exists)
    subtract
    write two files (total, and latest - previous)
    p.stat().st_mtime to find latest
    Args:
        results (): the results produced during this run
        singletons (bool): whether to read/write a singleton specific file

    Returns:
        the same results back-filtered to remove previous results
    """
    # try to pull out the historic results folder
    try:
        if (
            historic := get_config()['dataset_specific'].get('historic_results')
        ) is None:
            logging.info(
                'No `historic_results` key in config - no filtering took place'
            )
            return results
    except KeyError:
        logging.info('No `dataset_specific` key in config - no filtering took place')
        return results

    logging.info(f'filtering current results against data in {historic}')
    # get the latest result file from the folder
    latest_results = find_latest(
        results_folder=historic, start='singletons' if singletons else ''
    )

    logging.info(f'latest results: {latest_results}')

    if latest_results is None:
        # no results to subtract - current data IS cumulative data
        mini_results = make_cumulative_representation(results)

        save_new_cumulative(
            directory=historic, results=mini_results, singletons=singletons
        )

    else:
        cum_results = read_json_from_path(bucket_path=latest_results)

        # remove previously seen entries
        results = subtract_results(current=results, cumulative=cum_results)

        # add new entries into cumulative results
        add_results(results, cum_results)

        # save updated cumulative results
        save_new_cumulative(
            directory=historic, results=cum_results, singletons=singletons
        )

    return results


def make_cumulative_representation(results: dict[str, list[ReportedVariant]]) -> dict:
    """
    for the 'cumulative' representation, keep a minimised format
    the first time we generate a minimal format we need to translate

    Args:
        results (): the full results of the current run

    Returns:
        minimised form
    """

    mini_results = defaultdict(dict)
    for sample, variants in results.items():
        for var in variants:
            mini_results[sample][var.var_data.coords.string_format] = {
                'categories': {cat: TODAY for cat in var.var_data.categories},
                'support_vars': var.support_vars,
            }
    return mini_results


def save_new_cumulative(directory: str, results: dict, singletons: bool):
    """
    save the new cumulative results in the results dir
    include time & date in filename

    Args:
        directory ():
        results ():
        singletons (bool): whether to save this file as singleton-specific
    """

    new_file = to_path(directory) / f'{"singletons_" if singletons else ""}{TODAY}.json'
    with new_file.open('w') as handle:
        json.dump(results, handle, indent=4)
    logging.info(f'Wrote new cumulative data to {str(new_file)}')


def find_latest(results_folder: str, start: str = '', ext: str = 'json') -> str | None:
    """
    takes a directory of files, and finds the latest
    Args:
        results_folder (): local or remote folder
        start (str): the start of the filename, if applicable
        ext (): the type of files we're looking for

    Returns:
        timestamp-sorted latest file path, or None
    """

    date_sorted_files = sorted(
        to_path(results_folder).glob(f'{start}*.{ext}'),
        key=lambda x: x.stat().st_mtime,
        reverse=True,
    )
    if not date_sorted_files:
        return None

    return str(date_sorted_files[0].absolute())


def subtract_results(
    current: dict[str, list[ReportedVariant]], cumulative: dict
) -> dict[str, list[ReportedVariant]]:
    """
    take datasets of new and previous results (cumulative for this cohort)
    subtract all the previously seen results (exact)
        i.e. comp-het with a new partner should appear?
            - how does this gel with the category subtraction...
            - new partner = show all categories
    return the new-only

    Args:
        current (): results produced by this run
        cumulative (): results previously seen

    Returns:
        a diff between the two
    """

    # create a dict to contain novel results
    return_results = defaultdict(list)

    # iterate over all samples and their variants
    for sample, variants in current.items():

        # sample not seen before - take all variants
        if sample not in cumulative:
            return_results[sample] = variants
            continue

        # otherwise get ready to contain some variants
        return_results[sample] = []

        # check what has previously been reported for this sample
        for variant in variants:

            var_id = variant.var_data.coords.string_format

            # not seen - take in full
            if var_id not in cumulative[sample]:
                return_results[sample].append(variant)
                continue

            # if supported, allow if the support is new
            if variant.support_vars:

                # if novel supporting variants, don't need to check
                # categories, we want to show this. WALRUS
                if sups := [
                    sup
                    for sup in variant.support_vars
                    if sup not in cumulative[sample][var_id]['support_vars']
                ]:
                    variant.support_vars = sups
                    return_results[sample].append(variant)
                    continue

            # if seen, check for novel categories. Also, WALRUS
            if cats := [
                cat
                for cat in variant.var_data.categories
                if cat not in cumulative[sample][var_id]['categories'].keys()
            ]:
                # reduce the categories of interest for this round
                variant.var_data.categories = cats

                # and append the variant
                return_results[sample].append(variant)

    return dict(return_results)


def add_results(current: dict[str, list[ReportedVariant]], cumulative: dict):
    """
    take datasets of new and previous results (cumulative for this cohort)
    integrate the new results to form a new cumulative dataset
    the first time each variant-category is seen, add the current date

    Args:
        current (): results produced by this run
        cumulative (): results previously seen

    Returns:
        N/A - in place modification
        a union of all results, inc new categories for prev. variants
    """

    # iterate over all samples and their variants
    for sample, variants in current.items():

        # sample not seen before - take all variants
        if sample not in cumulative:
            cumulative[sample] = {}

        # each variant is an AbstractVariant
        for variant in variants:

            var_id = variant.var_data.coords.string_format

            # not seen - take in full
            if var_id not in cumulative[sample]:
                cumulative[sample][var_id] = {
                    'categories': {cat: TODAY for cat in variant.var_data.categories},
                    'support_vars': variant.support_vars,
                }

            else:

                # if seen, check for novel categories
                cumulative[sample][var_id]['categories'].update(
                    {
                        cat: TODAY
                        for cat in set(variant.var_data.categories)
                        - set(cumulative[sample][var_id]['categories'].keys())
                    }
                )
                cumulative[sample][var_id]['support_vars'] = sorted(
                    set(
                        variant.support_vars
                        + cumulative[sample][var_id]['support_vars']
                    )
                )
