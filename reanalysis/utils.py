"""
a collection of classes and methods
which may be shared across reanalysis components
"""


from typing import Any, Dict, List, Optional, Set, Tuple, Union
from dataclasses import dataclass, is_dataclass

import json
import logging
import re

import cyvcf2
from cloudpathlib import AnyPath
from cyvcf2 import Variant


# CompHetDict structure:
# {
#     sample: {
#         gene: {
#             variant: [variant, ...]
#         }
#     }
# }
# sample: string, e,g, CGP12345
# gene: string, e.g. ENSG012345
# variant: string, chr-pos-ref-alt
CompHetDict = Dict[str, Dict[str, Dict[str, List[str]]]]

InfoDict = Dict[str, Union[str, Dict[str, str]]]
PanelAppDict = Dict[str, Dict[str, Union[str, bool]]]

HOMREF: int = 0
HETALT: int = 1
UNKNOWN: int = 2
HOMALT: int = 3

# in cyVCF2, these ints represent HOMREF, and UNKNOWN
BAD_GENOTYPES: Set[int] = {HOMREF, UNKNOWN}


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
        var: Variant,
        samples: List[str],
        config: Dict[str, Any],
        as_singletons=False,
    ):

        # extract the coordinates into a separate object
        self.coords = Coordinates(
            var.CHROM.replace('chr', ''), var.POS, var.REF, var.ALT[0]
        )

        # set the class attributes
        self.category_1: bool = var.INFO.get('Category1') == 1
        self.category_2: bool = var.INFO.get('Category2') == 1
        self.category_3: bool = var.INFO.get('Category3') == 1

        # de novo category is a list of strings or empty list
        # if cohort runs as singletons, remove possibility of de novo
        if as_singletons:
            self.category_4 = []
        else:
            self.category_4: List[str] = (
                var.INFO.get('Category4').split(',')
                if var.INFO.get('Category4') != 'missing'
                else []
            )

        self.category_support: bool = var.INFO.get('CategorySupport') == 1

        # get all zygosities once per variant
        # abstraction avoids pulling per-sample calls again later
        self.het_samples, self.hom_samples = get_non_ref_samples(
            variant=var, samples=samples
        )

        self.info: Dict[str, str] = extract_info(variant=var)
        self.transcript_consequences: List[Dict[str, str]] = extract_csq(
            variant=var, config=config
        )

    @property
    def category_1_2_3(self) -> bool:
        """
        check that the variant has at least one assigned class
        supporting category is considered here
        :return:
        """
        return any(
            [
                self.category_1,
                self.category_2,
                self.category_3,
            ]
        )

    @property
    def category_non_support(self) -> bool:
        """
        check the variant has at least one non-support category assigned
        :return:
        """
        return any(
            [
                self.category_1,
                self.category_2,
                self.category_3,
                self.category_4,
            ]
        )

    @property
    def is_classified(self) -> bool:
        """
        check that the variant has at least one assigned class
        supporting category is considered here
        :return:
        """
        return any(
            [
                self.category_1,
                self.category_2,
                self.category_3,
                self.category_4,
                self.category_support,
            ]
        )

    @property
    def support_only(self) -> bool:
        """
        checks that the variant was only class 4
        :return:
        """
        return self.category_support and not any(
            [
                self.category_1,
                self.category_2,
                self.category_3,
                self.category_4,
            ]
        )

    def category_ints(self, sample: str) -> List[str]:
        """
        get a list of ints representing the classes present on this variant
        for each category, append that number if the class is present

        - support is not an int
        - de novo on a per-sample basis
        """
        categories = []

        # for the first 3 categories, append the value if flag is present
        for index, cat in enumerate(
            [self.category_1, self.category_2, self.category_3], 1
        ):
            if cat:
                categories.append(str(index))

        if self.sample_de_novo(sample_id=sample):
            categories.append('de_novo')

        if self.category_support:
            categories.append('in_silico')

        return categories

    def sample_de_novo(self, sample_id: str) -> bool:
        """
        takes a specific sample ID, to check if the sample has a de novo call

        :param sample_id:
        :return:
        """

        return sample_id in self.category_4

    def sample_specific_category_check(self, sample_id: str) -> bool:
        """

        :param sample_id:
        :return:
        """
        return self.category_1_2_3 or self.sample_de_novo(sample_id)


@dataclass
class ReportedVariant:
    """
    minimal model representing variant categorisation event
    the initial variant
    the MOI passed
    the support (if any)
    """

    sample: str
    gene: str
    var_data: AbstractVariant
    reasons: Set[str]
    supported: bool
    support_vars: Optional[List[str]] = None


def canonical_contigs_from_vcf(reader: cyvcf2.VCFReader) -> Set[str]:
    """
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
    variant_source: cyvcf2.VCFReader,
    config: Dict[str, Any],
    panelapp_data: PanelAppDict,
    singletons: bool,
    blacklist: Optional[List[str]] = None,
) -> Dict[str, Dict[str, AbstractVariant]]:
    """
    takes a cyvcf2.VCFReader instance, and a specified chromosome
    iterates over all variants in the region, and builds a lookup
    {
        gene: {
            var1_as_string: var1,
            var2_as_string: var2,
        },
        ...
    }
    :param contig: contig name from header (canonical_contigs_from_vcf)
    :param variant_source: VCF reader instance
    :param config: configuration file
    :param panelapp_data:
    :param singletons:
    :param blacklist:
    :return: populated lookup dict
    """

    if blacklist is None:
        blacklist = []

    # a dict to allow lookup of variants on this whole chromosome
    contig_variants = 0
    contig_dict = {}

    # iterate over all variants on this contig and store by unique key
    # if contig has no variants, prints an error and returns []
    for variant in variant_source(contig):
        abs_var = AbstractVariant(
            var=variant,
            samples=variant_source.samples,
            config=config,
            as_singletons=singletons,
        )

        if abs_var.coords.string_format in blacklist:
            logging.info(
                f'Skipping blacklisted variant: {abs_var.coords.string_format}'
            )
            continue

        # do category 2 'new' test
        # if the gene isn't 'new' in PanelApp, remove Class2 flag
        if abs_var.category_2:
            # implement the c2 check
            gene_id = abs_var.info.get('gene_id')
            gene_data = panelapp_data.get(gene_id, False)
            if not gene_data:
                continue

            # is the gene 'new'? if not, skip
            if not gene_data.get('new', False):
                variant.category_2 = False

        # if unclassified, skip the whole variant
        if not abs_var.is_classified:
            continue

        # update the variant count
        contig_variants += 1

        # update the gene index dictionary
        contig_dict.setdefault(abs_var.info.get('gene_id'), {})[
            abs_var.coords.string_format
        ] = abs_var

    logging.info(f'Contig {contig} contained {contig_variants} variants')
    logging.info(f'Contig {contig} contained {len(contig_dict)} genes')
    return contig_dict


def read_json_from_path(bucket_path: str) -> Dict[str, Any]:
    """
    take a path to a JSON file, read into an object
    :param bucket_path:
    """
    with open(AnyPath(bucket_path), encoding='utf-8') as handle:
        return json.load(handle)


def get_simple_moi(panel_app_moi: str) -> str:
    """
    takes the vast range of PanelApp MOIs, and reduces to a reduced
    range of cases which can be easily implemented in RD analysis
    This is required to reduce the complexity of an MVP
    Could become a strict enumeration

    {
        Biallelic: [all, panelapp, moi, equal, to, biallelic],
        Monoallelic: [all MOI to be interpreted as monoallelic]
    }

    :param panel_app_moi: full PanelApp string
    :return: a simplified representation
    """

    # default to considering both. NOTE! Many genes have Unknown MOI!
    simple_moi = 'Mono_And_Biallelic'
    if panel_app_moi is None or panel_app_moi == 'Unknown':
        # exit iteration, all simple moi considered
        return simple_moi

    # ideal for match-case, coming to a python 3.10 near you!
    if panel_app_moi.startswith('BIALLELIC'):
        simple_moi = 'Biallelic'
    if panel_app_moi.startswith('BOTH'):
        simple_moi = 'Mono_And_Biallelic'
    if panel_app_moi.startswith('MONO'):
        simple_moi = 'Monoallelic'
    if panel_app_moi.startswith('X-LINKED'):
        if 'biallelic' in panel_app_moi:
            simple_moi = 'Hemi_Bi_In_Female'
        else:
            simple_moi = 'Hemi_Mono_In_Female'

    return simple_moi


def get_non_ref_samples(
    variant: Variant, samples: List[str]
) -> Tuple[Set[str], Set[str]]:
    """
    for this variant, find all samples with a call
    cyvcf2 uses 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
    return het, hom, and the union of het and hom

    maybe something different would be more versatile
    e.g.
    hets = {
        sample: '0/1',
    }
    where the genotype and phased (|) vs unphased (/)
    is determined from the variant.genotypes attribute
    This would make it trivial for the final output to
    have an accurate representation of the parsed GT
    without having to regenerate the string rep.
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


def extract_csq(variant: Variant, config: Dict[str, Dict[str, str]]):
    """
    specifically handle extraction of the CSQ list
    :param variant:
    :param config:
    :return:
    """

    # config region concerning variant objects
    object_conf = config.get('variant_object')
    csq_string = object_conf.get('csq_string')

    # pull the CSQ content, and the format to decode it
    csq_contents = variant.INFO.get('CSQ')

    # break mono-CSQ-string into components
    csq_categories = list(map(str.lower, csq_string.split('|')))

    # iterate over all consequences, and make each into a dict
    return [
        dict(zip(csq_categories, each_csq.split('|')))
        for each_csq in csq_contents.split(',')
    ]


def extract_info(variant: Variant):
    """
    creates an INFO dict by pulling content from the variant info
    keeps a list of dictionaries for each transcript_consequence
    :param variant:
    :return:
    """

    # choose some values to exclude, and keep everything else
    exclusions = {
        'csq',
        'category1',
        'category2',
        'category3',
        'category4',
        'categorysupport',
        'support_only',
    }

    # grab the basic information from INFO
    return {x.lower(): y for x, y in variant.INFO if x.lower() not in exclusions}


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
