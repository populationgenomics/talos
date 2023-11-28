"""
A home for all data models used in AIP
"""
from enum import Enum

from pydantic import BaseModel, Field

from cpg_utils.config import get_config

from reanalysis.static_values import get_granular_date

NON_HOM_CHROM = ['X', 'Y', 'MT', 'M']
CHROM_ORDER = list(map(str, range(1, 23))) + NON_HOM_CHROM


class FileTypes(Enum):
    """
    enumeration of permitted input file types
    """

    HAIL_TABLE = '.ht'
    MATRIX_TABLE = '.mt'
    VCF = '.vcf'
    VCF_GZ = '.vcf.gz'
    VCF_BGZ = '.vcf.bgz'


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
        if self.chrom in CHROM_ORDER:
            return True
        return False

    def __eq__(self, other) -> bool:
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


class Variant(BaseModel):
    """
    the abstracted representation of a variant from any source
    """

    coordinates: Coordinates = Field(repr=True)
    info: dict[str, str | int | float | list[str] | bool] = Field(default_factory=dict)
    het_samples: set[str] = Field(default_factory=set, exclude=True)
    hom_samples: set[str] = Field(default_factory=set, exclude=True)
    boolean_categories: list[str] = Field(default_factory=list, exclude=True)
    sample_categories: list[str] = Field(default_factory=list, exclude=True)
    sample_support: list[str] = Field(default_factory=list, exclude=True)
    phased: dict = Field(default_factory=dict)

    def __str__(self):
        return repr(self)

    def __lt__(self, other):
        return self.coordinates < other.coordinates

    def __eq__(self, other):
        return self.coordinates == other.coordinates

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

    def sample_support_only(self, sample_id: str) -> bool:
        """
        check that the variant is exclusively cat. support
        check that this sample is missing from sample flags

        Returns:
            True if support only
        """
        return self.has_support and not (
            self.category_non_support or self.sample_categorised_check(sample_id)
        )

    def category_values(self, sample: str) -> list[str]:
        """
        get all variant categories
        steps category flags down to booleans - true for this sample

        Args:
            sample (str): sample id

        Returns:
            list of all categories applied to this variant
        """

        # step down all category flags to boolean flags
        categories = [
            category.replace('categorysample', '')
            for category in self.sample_categories
            if sample in self.info[category]  # type: ignore
        ]
        categories.extend(
            [
                bool_cat.replace('categoryboolean', '')
                for bool_cat in self.boolean_categories
                if self.info[bool_cat]
            ]
        )

        if self.has_support:
            categories.append('support')

        return categories

    def sample_categorised_check(self, sample_id: str) -> bool:
        """
        check if any *sample categories applied for this sample

        Args:
            sample_id (str): the specific sample ID to check

        Returns:
            bool: True if this sample features in any
                  named-sample category, includes 'all'
        """

        return any(
            sam in self.info[sam_cat]  # type: ignore
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
        big_cat = self.has_boolean_categories or self.sample_categorised_check(
            sample_id
        )
        if allow_support:
            return big_cat or self.has_support
        return big_cat


class SmallVariant(Variant):
    depths: dict[str, int] = Field(default_factory=dict, exclude=True)
    ab_ratios: dict[str, float] = Field(default_factory=dict, exclude=True)
    transcript_consequences: list[dict[str, str | float | int]] = Field(
        default_factory=list
    )

    def get_sample_flags(self, sample: str) -> list[str]:
        """
        gets all report flags for this sample - currently only one flag
        """
        return self.check_ab_ratio(sample) + self.check_read_depth(sample)

    def check_read_depth(self, sample: str, threshold: int = 10) -> list[str]:
        """
        flag low read depth for this sample

        Args:
            sample (str): sample ID matching VCF
            threshold (int): cut-off for flagging as a failure

        Returns:
            return a flag if this sample has low read depth
        """
        if self.depths[sample] < threshold:
            return ['Low Read Depth']
        return []

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


class StructuralVariant(Variant):
    def check_ab_ratio(self, *args, **kwargs) -> list[str]:
        """
        dummy method for AB ratio checking - not implemented for SVs
        """
        return []

    def get_sample_flags(self, *args, **kwargs) -> list[str]:
        """
        dummy method for flag checking - not implemented for SVs (yet)
        """
        return []

    def check_read_depth(self, *args, **kwargs) -> list[str]:
        """
        dummy method for read depth checking - not implemented for SVs
        """
        return []


class ReportPanel(BaseModel):
    """
    simple storage for all the panels to present in tooltips
    """

    forced: list[str] = Field(default_factory=list)
    matched: list[str] = Field(default_factory=list)


class ReportVariant(BaseModel):
    """
    A variant passing MOI tests, to be reported
    """

    sample: str
    var_data: SmallVariant | StructuralVariant
    categories: list[str] = Field(default_factory=list)
    family: str = Field(default_factory=str)
    first_seen: str = Field(default=get_granular_date())
    flags: list[str] = Field(default_factory=list)
    gene: str = Field(default_factory=str)
    genotypes: dict[str, str] = Field(default_factory=dict)
    independent: bool = Field(default=False)
    labels: list[str] = Field(default_factory=list)
    panels: ReportPanel = Field(default_factory=ReportPanel)
    phenotypes: list[str] = Field(default_factory=list)
    reasons: set[str] = Field(default_factory=set)
    support_vars: set[str] = Field(default_factory=set)

    def __eq__(self, other):
        """
        makes reported variants comparable
        """
        # self_supvar = set(self.support_vars)
        # other_supvar = set(other.support_vars)
        return (
            self.sample == other.sample
            and self.var_data.coordinates == other.var_data.coordinates
        )

    def __lt__(self, other):
        return self.var_data.coordinates < other.var_data.coordinates


class PanelDetail(BaseModel):
    """
    A gene from PanelApp, combining all MOI and panel IDs
    where the gene features on multiple panels
    """

    symbol: str
    chrom: str = Field(default_factory=str)
    all_moi: set[str] = Field(default_factory=set)
    moi: str = Field(default_factory=str)
    new: list[int] = Field(default_factory=list)
    panels: list[int] = Field(default_factory=list)


class PanelShort(BaseModel):
    """
    Short panel summary, used in the metadata section
    """

    name: str
    version: str
    id: int


class PanelList(BaseModel):
    panels: list[PanelShort] = Field(default_factory=list)


class PanelApp(BaseModel):
    metadata: list[PanelShort]
    genes: dict[str, PanelDetail]


class HistoricPanels(BaseModel):
    ...


class CategoryMeta(BaseModel):
    """
    The mapping of category names to their display names
    """

    categories: dict[str, str] = Field(default=get_config()['categories'])


class HistoricSampleVariant(BaseModel):
    """ """

    categories: dict[str, str]
    support_vars: list[str] = Field(
        default_factory=list,
        description='supporting variants if this has been identified in a comp-het',
    )
    independent: bool = Field(default=True)


class HistoricVariants(BaseModel):
    """
    The model representing the state transition file
    All relevant metadata relating to the available categories
    Then a per-participant dict of variants, containing the categories
    they have been assigned, date first seen, and supporting variants
    """

    metadata: CategoryMeta = Field(default_factory=CategoryMeta)
    results: dict[str, HistoricSampleVariant] = Field(default_factory=dict)


class ResultMeta(BaseModel):
    """
    metadata for a result set
    """

    input_file: str
    cohort: str
    run_datetime: str = Field(default=get_granular_date())
    family_breakdown: dict[str, int] = Field(default_factory=dict)
    panels: PanelList = Field(default_factory=PanelList)
    container: str = Field(default_factory=str)
    categories: dict[str, str] = Field(default=get_config()['categories'])


class MemberSex(Enum):
    MALE = 'male'
    FEMALE = 'female'
    UNKNOWN = 'unknown'


class FamilyMembers(BaseModel):
    affected: bool = Field(default=False)
    ext_id: str = Field(default_factory=str)
    sex: MemberSex = Field(default=MemberSex.UNKNOWN)


class ParticipantMeta(BaseModel):
    ext_id: str
    family_id: str
    members: dict[str, FamilyMembers] = Field(default_factory=dict)
    phenotypes: list[str] = Field(default_factory=list)
    panel_ids: list[int] = Field(default_factory=list)
    panel_names: list[str] = Field(default_factory=list)
    solved: bool = Field(default=False)


class ParticipantResults(BaseModel):
    """
    A representation of a result set
    """

    variants: list[ReportVariant]
    metadata: ParticipantMeta


class ResultData(BaseModel):
    """
    A representation of a result set
    """

    results: dict[str, ParticipantResults]
    metadata: ResultMeta


class ModelVariant(BaseModel):
    """
    might be required for the VCF generator
    """

    ...


class PhenoMatchedPanels(BaseModel):
    ...
