"""
a collection of classes and methods
which may be shared across reanalysis components
"""


from typing import Any, Dict, List, Optional, Set, Tuple, Union
from dataclasses import dataclass, is_dataclass
from csv import DictReader
import json
import re
import sys

import cyvcf2
from cloudpathlib import AnyPath
from cyvcf2 import Variant

# modify how to import cache based on the system python
v_info = sys.version_info
if v_info.major != 3:
    raise Exception('Come on, use Python3')
if v_info.minor < 9:
    from functools import lru_cache as cache
else:
    from functools import cache

# sample: gene: variant: [partner_variant, ...]
CompHetDict = Dict[str, Dict[str, Dict[str, List[str]]]]
InfoDict = Dict[str, Union[str, Dict[str, str]]]
PanelAppDict = Dict[str, Dict[str, Union[str, bool]]]

COMP_HET_VALUES: List[str] = ['sample', 'gene', 'id', 'chrom', 'pos', 'ref', 'alt']
VARIANT_STRING_TEMPLATE: str = '{}-{}-{}-{}'
COMP_HET_TEMPLATE: str = f'{VARIANT_STRING_TEMPLATE}-{{}}'

HOMREF: int = 0
HETALT: int = 1
UNKNOWN: int = 2
HOMALT: int = 3

# in cyVCF2, these ints represent HOMREF, and UNKNOWN
BAD_GENOTYPES: Set[int] = {HOMREF, UNKNOWN}


def canonical_contigs_from_vcf(reader: cyvcf2.VCFReader) -> Set[str]:
    """
    read the header fields from the VCF handle
    return a set of all 'canonical' contigs
    :param reader:
    :return:
    """

    # contig matching regex - remove all HLA/decoy/unknown
    contig_re = re.compile(r'^(chr)?[1-9XYMT]{1,2}$')

    return {
        contig['ID']
        for contig in reader.header_iter()
        if contig['HeaderType'] == 'CONTIG' and re.match(contig_re, contig['ID'])
    }


def read_json_dict_from_path(bucket_path: str) -> Dict[str, Any]:
    """
    take a path to a JSON file, read into an object
    :param bucket_path:
    """
    with open(AnyPath(bucket_path), encoding='utf-8') as handle:
        return json.load(handle)


@cache
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


@dataclass
class PedPerson:
    """
    holds attributes about a single PED file entry
    this will need to be enhanced for family analysis
    a prototype already exists in the prototype folder
    """

    sample: str
    male: bool
    affected: bool


def parse_ped_simple(ped: str) -> Dict[str, PedPerson]:
    """
    take individual attributes - sample ID, sex, affected
    :param ped: path to the ped file
    :return:
    """

    ped_dict: Dict[str, PedPerson] = {}
    with open(AnyPath(ped), 'r', encoding='utf-8') as handle:
        for line in DictReader(handle, delimiter='\t'):

            # slot in the sample ID and two Booleans
            ped_dict[line['Individual ID']] = PedPerson(
                line['Individual ID'], line['Sex'] == '1', line['Affected'] == '2'
            )
    return ped_dict


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


def extract_info(variant: Variant, config: InfoDict):
    """
    creates an INFO dict by pulling content from the variant info
    keeps a list of dictionaries for each transcript_consequence
    :param variant:
    :param config:
    :return:
    """

    # choose some values to exclude, and keep everything else
    exclusions = {'CSQ', 'Class1', 'Class2', 'Class3', 'Class4'}

    # config region concerning variant objects
    object_conf = config.get('variant_object')

    # grab the basic information from INFO
    info_dict = {x.lower(): y for x, y in variant.INFO if x not in exclusions}

    # pull the CSQ content, and the format to decode it
    csq_contents = variant.INFO.get('CSQ')
    csq_string = object_conf.get('csq_string')

    # break mono-CSQ-string into components
    csq_categories = list(map(str.lower, csq_string.split('|')))

    # iterate over all consequences, and make each into a dict
    info_dict['transcript_consequences'] = [
        dict(zip(csq_categories, each_csq.split('|')))
        for each_csq in csq_contents.split(',')
    ]
    return info_dict


@dataclass
class Coordinates:
    """
    a home for the extracted variant attributes
    """

    chrom: str
    pos: int
    ref: str
    alt: str

    @property
    def string_format(self):
        """
        forms a string representation
        """
        return VARIANT_STRING_TEMPLATE.format(self.chrom, self.pos, self.ref, self.alt)


@dataclass
class AbstractVariant:  # pylint: disable=too-many-instance-attributes
    """
    create a bespoke variant class
    pull all content out of the cyvcf2 object
    we could have a separate implementation for pyvcf,
    or for a direct parser...
    """

    def __init__(self, var: Variant, samples: List[str], config: Dict[str, Any]):

        # extract the coordinates into a separate object
        self.coords = Coordinates(
            var.CHROM.replace('chr', ''), var.POS, var.REF, var.ALT[0]
        )

        # set the class attributes
        self.category_1: bool = var.INFO.get('Category1') == 1
        self.category_2: bool = var.INFO.get('Category2') == 1
        self.category_3: bool = var.INFO.get('Category3') == 1
        self.category_4: bool = var.INFO.get('Category4') == 1

        # get all zygosities once per variant
        # abstraction avoids pulling per-sample calls again later
        self.het_samples, self.hom_samples = get_non_ref_samples(
            variant=var, samples=samples
        )

        self.info: InfoDict = extract_info(variant=var, config=config)

    @property
    def is_classified(self) -> bool:
        """
        check that the variant has at least one assigned class
        :return:
        """
        return any([self.category_1, self.category_2, self.category_3, self.category_4])

    @property
    def category_4_only(self) -> bool:
        """
        checks that the variant was only class 4
        :return:
        """
        return self.category_4 and not any(
            [
                self.category_1,
                self.category_2,
                self.category_3,
            ]
        )

    @property
    def category_ints(self) -> List[int]:
        """
        get a list of ints representing the classes present on this variant
        for each numerical class, append that number if the class is present
        """
        return [
            integer
            for integer, category_bool in enumerate(
                [
                    self.category_1,
                    self.category_2,
                    self.category_3,
                    self.category_4,
                ],
                1,
            )
            if category_bool
        ]


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
