"""
A home for all data models used in Talos
"""

from enum import Enum
from itertools import pairwise
from typing import Any

from loguru import logger
from pydantic import BaseModel, Field

from talos.liftover.lift_1_0_0_to_1_0_1 import resultdata as rd_100_to_101
from talos.liftover.lift_1_0_2_to_1_0_3 import resultdata as rd_102_to_103
from talos.liftover.lift_1_0_3_to_1_1_0 import resultdata as rd_103_to_110
from talos.liftover.lift_1_1_0_to_1_2_0 import resultdata as rd_110_to_120
from talos.liftover.lift_1_2_0_to_2_0_0 import panelapp as pa_120_to_200
from talos.liftover.lift_1_2_0_to_2_0_0 import resultdata as rd_120_to_200
from talos.liftover.lift_2_0_0_to_2_1_0 import panelapp as pa_200_to_210
from talos.liftover.lift_2_0_0_to_2_1_0 import resultdata as rd_200_to_210
from talos.liftover.lift_2_1_0_to_2_2_0 import dl_panelapp as dl_pa_210_to_220
from talos.liftover.lift_2_1_0_to_2_2_0 import resultdata as rd_210_to_220
from talos.liftover.lift_none_to_1_0_0 import resultdata as rd_none_to_1_0_0
from talos.static_values import get_granular_date

NON_HOM_CHROM = ['X', 'Y', 'MT', 'M']
CHROM_ORDER = list(map(str, range(1, 23))) + NON_HOM_CHROM

# some kind of version tracking
CURRENT_VERSION = '2.2.0'
ALL_VERSIONS = [None, '1.0.0', '1.0.1', '1.0.2', '1.0.3', '1.1.0', '1.2.0', '2.0.0', '2.1.0', '2.2.0']

# ratios for use in AB testing
MAX_WT = 0.15
MIN_HET = 0.25
MAX_HET = 0.75
MIN_HOM = 0.85

# dictionary for translating all the previous aliases for categories to their more descriptive names
CATEGORY_TRANSLATOR: dict[str, str] = {
    '1': 'ClinVar P/LP',
    'clinvarplp': 'ClinVar P/LP',
    'ClinVarP/LP': 'ClinVar P/LP',
    'clinvar0star': 'ClinVar 0-star',
    'clinvar0starnewgene': 'ClinVar Recent Gene',
    '3': 'High Impact',
    'highimpact': 'High Impact',
    '4': 'De Novo',
    'denovo': 'De Novo',
    '5': 'SpliceAI',
    'spliceai': 'SpliceAI',
    '6': 'AlphaMissense',
    'alphamissense': 'AlphaMissense',
    'pm5': 'PM5',
    'sv1': 'LOF SV',
    'lofsv': 'LOF SV',
    'svdb': 'SpliceVarDB',
    'splicevardb': 'SpliceVarDB',
    'exomiser': 'Exomiser',
}


def translate_category(cat: str) -> str:
    """Translate a category from config file to a more descriptive name. If not found, return the original."""
    return CATEGORY_TRANSLATOR.get(cat.lower(), cat)


class FileTypes(Enum):
    """
    enumeration of permitted input file types
    """

    HAIL_TABLE = '.ht'
    MATRIX_TABLE = '.mt'
    PED = 'ped'
    VCF = '.vcf'
    VCF_GZ = '.vcf.gz'
    VCF_BGZ = '.vcf.bgz'


class HpoTerm(BaseModel):
    """
    A representation of a HPO term
    """

    id: str
    label: str


class Coordinates(BaseModel):
    """
    A representation of genomic coordinates
    """

    chrom: str
    pos: int
    ref: str
    alt: str

    @property
    def string_format(self) -> str:
        """
        forms a string representation: chr-pos-ref-alt
        """
        return f'{self.chrom}-{self.pos}-{self.ref}-{self.alt}'

    def __lt__(self, other) -> bool:
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
        return self.chrom in CHROM_ORDER


class VariantCommon(BaseModel):
    """
    the abstracted representation of a variant from any source
    """

    coordinates: Coordinates = Field(repr=True)
    info: dict[str, Any] = Field(default_factory=dict)
    het_samples: set[str] = Field(default_factory=set, exclude=True)
    hom_samples: set[str] = Field(default_factory=set, exclude=True)
    boolean_categories: list[str] = Field(default_factory=list, exclude=True)
    sample_categories: list[str] = Field(default_factory=list, exclude=True)
    ignored_categories: set[str] = Field(default_factory=set)
    support_categories: set[str] = Field(default_factory=set)
    phased: dict = Field(default_factory=dict, exclude=True)

    def __str__(self):
        return repr(self)

    def __lt__(self, other):
        return self.coordinates < other.coordinates

    def __eq__(self, other):
        return self.coordinates == other.coordinates

    def __hash__(self):
        return hash(self.coordinates)

    def category_values(self, sample: str) -> set[str]:
        """
        get all variant categories which apply to this Sample (Boolean and Sample-specific)
        steps category flags down to booleans - true for this sample

        Args:
            sample (str): sample id

        Returns:
            set of all categories applied to this variant, using an enhanced labelling String for reporting purposes.
        """

        # step down all category flags to boolean flags
        categories: set[str] = set()
        for category in self.sample_categories:
            cat_samples = self.info[category]
            if not isinstance(cat_samples, list):
                raise TypeError(f'Sample categories should be a list: {cat_samples}')

            if any(sam in cat_samples for sam in [sample, 'all']):
                categories.add(category.removeprefix('categorysample'))

        categories.update(
            {bool_cat.replace('categoryboolean', '') for bool_cat in self.boolean_categories if self.info[bool_cat]},
        )

        # upgrade all the category labels for the report
        return {translate_category(cat) for cat in categories}

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

        # get all categories, both boolean and sample-specific - set of Strings
        categories_applied = self.category_values(sample_id)

        # if we don't want to allow support variants, remove any from the list of applied categories
        if not allow_support:
            # add the longer names to the support_categories - this is a workaround for the fact that the support
            # categories entry in the config file can now be the numberical/short IDs, or longer names
            remove_support = self.support_categories
            remove_support.update({translate_category(x) for x in self.support_categories})
            categories_applied -= remove_support

        return len(categories_applied) > 0

    def dodgy_ab_ratio_test(self, *args, **kwargs) -> set[str]:  # noqa: ARG002, ANN002, ANN003
        """
        dummy method for AB ratio checking - not implemented for SVs
        """
        return set()

    def get_sample_flags(self, *args, **kwargs) -> set[str]:  # noqa: ARG002, ANN002, ANN003
        """
        dummy method for flag checking - not implemented for SVs (yet)
        """
        return set()

    def insufficient_read_depth(self, *args, **kwargs) -> set[str]:  # noqa: ARG002, ANN002, ANN003
        """
        dummy method for read depth checking - not implemented for SVs
        """
        return set()

    def insufficient_alt_depth(self, sample: str, threshold: int = 5) -> set[str]:  # noqa: ARG002
        """
        dummy method for alt read depth checking - not implemented for SVs
        """
        return set()

    def min_alt_ratio(self, sample: str, threshold: float = 0.1) -> bool:  # noqa: ARG002
        """Dummy method for alt read ratio checking - not implemented for SVs."""
        return True


class SmallVariant(VariantCommon):
    """
    a representation of a small variant
    note that transcript_consequences is not optional
    we require that something specific to SmallVariant(s) is mandatory
    this is in order to correctly deserialise the Small/Structural objects
    into the appropriate types. If it is optional, Pydantic can coerce
    everything as a SmallVariant
    """

    alt_depths: dict[str, int] = Field(default_factory=dict, exclude=True)
    depths: dict[str, int] = Field(default_factory=dict, exclude=True)
    ab_ratios: dict[str, float] = Field(default_factory=dict, exclude=True)
    transcript_consequences: list[dict[str, str | float | int]]

    def get_sample_flags(self, sample: str) -> set[str]:
        """
        gets all report flags for this sample - currently only one flag
        """
        return self.dodgy_ab_ratio_test(sample) | self.insufficient_read_depth_strict(sample)

    def insufficient_read_depth(self, sample: str, threshold: int = 10, var_is_cat_1: bool = False) -> set[str]:
        """
        flag low read depth for this sample
        this version is used to determine whether a variant should be considered for _this_ sample
        in this situation we pass variants in clinvar, regardless of read depth

        Args:
            sample (str): sample ID matching VCF
            threshold (int): cut-off for flagging as a failure
            var_is_cat_1 (bool): flag if this variant is a category 1

        Returns:
            return a flag if this sample has low read depth
        """
        if var_is_cat_1:
            return set()
        return self.insufficient_read_depth_strict(sample, threshold)

    def insufficient_read_depth_strict(self, sample: str, threshold: int = 10) -> set[str]:
        """
        flag low read depth for this sample - doesn't care if it's in ClinVar

        Args:
            sample (str): sample ID matching VCF
            threshold (int): cut-off for flagging as a failure

        Returns:
            return a flag if this sample has low read depth
        """
        if self.depths[sample] < threshold:
            return {'Low Read Depth'}
        return set()

    def insufficient_alt_depth(self, sample: str, threshold: int = 5):
        """
        flag variants which have insufficient alt read support

        Args:
            sample (str): sample ID matching the VCF
            threshold (int): minimum number of alt reads required

        Returns:
            return a flag if this sample has low supporting read depth
        """
        if self.alt_depths[sample] < threshold:
            return {'Low Alt Support'}
        return set()

    def dodgy_ab_ratio_test(self, sample: str) -> set[str]:
        """
        AB ratio test for this sample's variant call

        Args:
            sample (str): sample ID

        Returns:
            set[str]: empty, or indicating an AB ratio failure
        """
        het = sample in self.het_samples
        hom = sample in self.hom_samples
        variant_ab = self.ab_ratios.get(sample, 0.0)

        if (variant_ab <= MAX_WT) or (het and not MIN_HET <= variant_ab <= MAX_HET) or (hom and variant_ab <= MIN_HOM):
            return {'AB Ratio'}
        return set()

    def min_alt_ratio(self, sample: str, threshold: float = 0.2) -> bool:
        sample_depth = self.depths[sample]
        sample_alt = self.alt_depths[sample]

        return (sample_alt / sample_depth) > threshold


class StructuralVariant(VariantCommon):
    """
    placeholder for any methods/data specific to Structural Variants
    """


# register all interchangeable models here
VARIANT_MODELS = SmallVariant | StructuralVariant


class ReportPanel(BaseModel):
    """
    simple storage for all the panels to present in tooltips
    """

    forced: dict[int, str] = Field(default_factory=dict)
    matched: dict[int, str] = Field(default_factory=dict)


class ReportVariant(BaseModel):
    """
    A variant passing MOI tests, to be reported
    """

    sample: str
    var_data: VARIANT_MODELS
    categories: dict[str, str] = Field(default_factory=dict)
    date_of_phenotype_match: str | None = None

    phenotype_labels: set[str] = Field(default_factory=set)

    evidence_last_updated: str = Field(default=get_granular_date())

    family: str = Field(default_factory=str)
    # 'tagged' is seqr-compliant language
    first_tagged: str = Field(default=get_granular_date())
    flags: set[str] = Field(default_factory=set)
    gene: str = Field(default_factory=str)
    genotypes: dict[str, str] = Field(default_factory=dict)
    labels: set[str] = Field(default_factory=set)
    panels: ReportPanel = Field(default_factory=ReportPanel)
    phenotypes: list[HpoTerm] = Field(default_factory=list)
    reasons: str = Field(default_factory=str)
    support_vars: set[str] = Field(default_factory=set)

    # new, recording this here instead of in the history file
    clinvar_stars: int | None = None

    # log whether there was an increase in ClinVar star rating since the last run
    clinvar_increase: bool = Field(default=False)

    # exomiser results - I'd like to store this in a cleaner way in future
    exomiser_results: list[str] = Field(default_factory=list)
    found_in_current_run: bool = Field(default=True)

    def __eq__(self, other):
        """
        makes reported variants comparable
        """
        return self.sample == other.sample and self.var_data.coordinates == other.var_data.coordinates

    def __hash__(self):
        return hash((self.sample, self.var_data.coordinates))

    def __lt__(self, other):
        return self.var_data.coordinates < other.var_data.coordinates


class ParticipantHPOPanels(BaseModel):
    family_id: str = Field(default_factory=str)
    hpo_terms: list[HpoTerm] = Field(default_factory=list)
    panels: set[int] = Field(default_factory=set)
    matched_genes: set[str] = Field(default_factory=set)
    matched_phenotypes: set[str] = Field(default_factory=set)


class PanelDetail(BaseModel):
    """
    A gene from PanelApp, combining all MOI and panel IDs
    where the gene features on multiple panels
    """

    symbol: str
    chrom: str = Field(default_factory=str)
    moi: str = Field(default_factory=str)
    new: set[int] = Field(default_factory=set)
    panels: set[int] = Field(default_factory=set)


class PanelShort(BaseModel):
    """
    Short panel summary, used in the metadata section
    This object contains the coarse panel details
    """

    id: int
    name: str = Field(default_factory=str)
    version: str = 'UNKNOWN'


class PanelApp(BaseModel):
    # the PanelShort object contains id, but we use this in a few places to search for the name/version of a panel by id
    # having this as a dictionary of {id: {id: X, name: Y}} looks a bit wasteful, but simplifies code in a few places
    metadata: dict[int, PanelShort] = Field(default_factory=dict)
    genes: dict[str, PanelDetail] = Field(default_factory=dict)
    participants: dict[str, ParticipantHPOPanels] = Field(default_factory=dict)
    version: str = CURRENT_VERSION
    creation_date: str = Field(default=get_granular_date())


class DownloadedPanelAppGenePanelDetail(BaseModel):
    """ """

    moi: str
    date: str = Field(default_factory=str)


class DownloadedPanelAppGene(BaseModel):
    """ """

    symbol: str
    chrom: str = Field(default_factory=str)
    mane_symbol: str = Field(default_factory=str)
    ensg: str = Field(default_factory=str)
    # for every panel this gene has featured in, when did it become Green, and what was the MOI
    panels: dict[int, DownloadedPanelAppGenePanelDetail] = Field(default_factory=dict)


class DownloadedPanelApp(BaseModel):
    """ """

    # all panels and versions
    versions: list[PanelShort] = Field(default_factory=list)
    genes: dict[str, DownloadedPanelAppGene] = Field(default_factory=dict)
    hpos: dict[int, list[HpoTerm]] = Field(default_factory=dict)
    version: str = CURRENT_VERSION
    date: str = Field(default=get_granular_date())


class ResultMeta(BaseModel):
    """
    metadata for a result set
    """

    version: str = Field(default_factory=str)
    family_breakdown: dict[str, int] = Field(default_factory=dict)
    input_file: str = Field(default_factory=str)
    panels: dict[int, PanelShort] = Field(default_factory=dict)
    run_datetime: str = Field(default=get_granular_date())

    # a count of variants per category, used for the report
    variant_breakdown: dict[str, dict[str, float | int]] = Field(default_factory=dict)
    # list of samples that had no variants in this result set
    samples_with_no_variants: list[str] = Field(default_factory=list)
    # external labels provided but not matched to any variant in this result set
    unused_ext_labels: list[dict[str, Any]] = Field(default_factory=list)


class MemberSex(Enum):
    UNKNOWN = 'unknown'
    MALE = 'male'
    FEMALE = 'female'


class FamilyMembers(BaseModel):
    affected: bool = Field(default=False)
    sex: str = Field(default=MemberSex.UNKNOWN.value)


class ParticipantMeta(BaseModel):
    family_id: str
    members: dict[str, FamilyMembers] = Field(default_factory=dict)
    phenotypes: list[HpoTerm] = Field(default_factory=list)
    panel_details: dict[int, PanelShort] = Field(default_factory=dict)
    solved: bool = Field(default=False)
    present_in_small: bool = Field(default=False)
    present_in_sv: bool = Field(default=False)


class ParticipantResults(BaseModel):
    """
    A representation of a result set
    """

    variants: list[ReportVariant] = Field(default_factory=list)
    metadata: ParticipantMeta = Field(default_factory=ParticipantMeta)


class ResultData(BaseModel):
    """
    A representation of a result set
    """

    results: dict[str, ParticipantResults] = Field(default_factory=dict)
    metadata: ResultMeta = Field(default_factory=ResultMeta)
    version: str = CURRENT_VERSION


class MiniVariant(BaseModel):
    categories: set[str] = Field(default_factory=set)
    support_vars: set[str] = Field(default_factory=set)
    independent: bool = Field(default=True)


class MiniForSeqr(BaseModel):
    results: dict[str, dict[str, MiniVariant]] = Field(default_factory=dict)


class PedigreeMember(BaseModel):
    """
    This will be a more searchable implementation of the peds pedigree
    """

    family: str
    id: str
    mother: str | None = None
    father: str | None = None
    sex: int
    affected: int
    ext_id: str = 'Missing'
    hpo_terms: set[str] = Field(default_factory=set)

    def __str__(self):
        """String representation of the participant, used when writing a Pedigree."""
        return f'{self.family}\t{self.id}\t{self.father}\t{self.mother}\t{self.sex}\t{self.affected}'


# methods defining how to transition between model versions. If unspecified, no transition is required
LIFTOVER_METHODS: dict = {
    DownloadedPanelApp: {
        '2.1.0_2.2.0': dl_pa_210_to_220,
    },
    PanelApp: {
        '1.2.0_2.0.0': pa_120_to_200,
        '2.0.0_2.1.0': pa_200_to_210,
    },
    ResultData: {
        'None_1.0.0': rd_none_to_1_0_0,
        '1.0.0_1.0.1': rd_100_to_101,
        '1.0.2_1.0.3': rd_102_to_103,
        '1.0.3_1.1.0': rd_103_to_110,
        '1.1.0_1.2.0': rd_110_to_120,
        '1.2.0_2.0.0': rd_120_to_200,
        '2.0.0_2.1.0': rd_200_to_210,
        '2.1.0_2.2.0': rd_210_to_220,
    },
}


def lift_up_model_version(
    data: dict,
    model: ResultData | PanelApp | DownloadedPanelApp,
) -> dict:
    """
    lift over data from one version to another
    this takes a dictionary object, and the model it is supposed to be
    we walk up from the version in the dict to the current models.py one,
    and apply all required transitions. This can be type changes, removing/adding mandatory fields, etc.
    incremental transition methods should be registered in reanalysis/liftover/old_to_new.py

    We expect that the majority of transitions will be no-ops,
    but we need to check for the existence of a method

    Args:
        data (dict): the model data prior to any transitions
        model (class): the data model the dict needs to be parsed as

    Returns:
        the input dictionary, transitioned to current format
    """

    # get current version, can be None (initial models did not specify)
    from_version = data.get('version')

    if from_version == CURRENT_VERSION:
        return data

    if from_version not in ALL_VERSIONS:
        raise ValueError(f'Unknown {model.__name__} version: {from_version}')

    if from_version != CURRENT_VERSION and not LIFTOVER_METHODS.get(model):
        logger.info(f'No liftover methods for {model.__name__}')

    # now pairwise iterate over ALL_VERSIONS, and find liftovers between each
    from_version_index = ALL_VERSIONS.index(from_version)

    # liftover pairs are (version, next_version)
    # e.g. [(None, 1), (1,2), (2,3), (3,4), (4,5)]
    # starting at the current index, and moving up from there
    # e.g. liftover from 2 to 5 would use [(2,3), (3,4), (4,5)]
    for previous, current in list(pairwise(ALL_VERSIONS))[from_version_index:]:
        liftover_key = f'{previous}_{current}'
        if liftover_key in LIFTOVER_METHODS[model]:
            data = LIFTOVER_METHODS[model][liftover_key](data)
        data['version'] = current
    return data
