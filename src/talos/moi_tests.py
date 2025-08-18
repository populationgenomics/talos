"""
A number of classes, each representing one Mode of Inheritance
One class (MoiRunner) to run all the appropriate MOIs on a variant

Reduce the PanelApp plain text MOI description into a few categories
We then run a permissive MOI match for the variant

Expected available gnomad annotations:
'gnomad': struct {
    gnomad_AC: int32,
    gnomad_AF: float64,
    gnomad_AC_XY: int32,
    gnomad_HomAlt: int32
}
"""

# mypy: ignore-errors
from abc import abstractmethod
from dataclasses import dataclass
from typing import ClassVar

from talos.config import config_retrieve
from talos.models import VARIANT_MODELS, ReportVariant, SmallVariant, StructuralVariant
from talos.pedigree_parser import PedigreeParser
from talos.utils import X_CHROMOSOME, CompHetDict

HEMI_CHROMS = {'chrX, chrY'}
SV_HEMI = {'male_n_hemialt'}
SV_HOMS = {'male_n_homalt', 'female_n_homalt'}


@dataclass
class GlobalFilter:
    """
    A Filter class, used to apply to any non-ClinVar Pathogenic Variants
    This pulls in a number of thresholds from the config file, and contains a 'too_common' method
    A variant run through this method will return true based on the thresholds if it's 'too_common'
    """

    # minimum variant AC to run callset frequency filters
    ac_threshold: ClassVar[int] = config_retrieve(['ValidateMOI', 'min_callset_ac_to_filter'])
    small_af: ClassVar[float] = config_retrieve(['ValidateMOI', 'callset_max_af'])
    sv_af: ClassVar[float] = config_retrieve(['ValidateMOI', 'callset_sv_max_af'])

    # a lookup of the attribute name vs. the corresponding configurable filter to be used on small variants
    small_dict: ClassVar[dict[str, float | int]] = {
        'gnomad_af': config_retrieve(['ValidateMOI', 'gnomad_max_af']),
        'gnomad_homalt': config_retrieve(['ValidateMOI', 'gnomad_max_homozygotes']),
    }

    # only to be applied on chrX/Y
    small_gnomad_hemi: ClassVar[int] = config_retrieve(['ValidateMOI', 'gnomad_max_hemizygotes'])

    # filters specific to SVs
    sv_dict: ClassVar[dict[str, float]] = {
        'gnomad_v2.1_sv_AF': config_retrieve(['ValidateMOI', 'gnomad_sv_max_af']),
    }

    def too_common(self, variant: SmallVariant | StructuralVariant) -> bool:
        """
        Check if a variant is too common in the population

        Args:
            variant (SmallVariant | StructuralVariant): the variant to check

        Returns:
            bool: True if the variant is too common
        """

        # check against each small-variant filter
        if isinstance(variant, SmallVariant):
            for key, threshold in self.small_dict.items():
                if key in variant.info and variant.info[key] > threshold:
                    return True

            # if there are sufficient instances, check for frequency in the callset
            if variant.info['ac'] > self.ac_threshold and variant.info['af'] > self.small_af:
                return True

            # on sex chroms, apply hemi-count filter
            if variant.coordinates.chrom in HEMI_CHROMS:
                return variant.info.get('gnomad_ac_xy', 0) > self.small_gnomad_hemi

        # check against the SV filters
        elif isinstance(variant, StructuralVariant):
            for key, threshold in self.sv_dict.items():
                if key in variant.info and variant.info[key] > threshold:
                    return True

            # if there are sufficient instances, check for frequency in the callset
            if variant.info['ac'] > self.ac_threshold and variant.info['af'] > self.sv_af:
                return True

        else:
            raise ValueError('Variant type not recognised')

        return False


@dataclass
class DominantFilter:
    """
    Similar to the GlobalFilter, but with stricter thresholds
    This is designed to run on variants being considered for Dominant inheritance
    """

    # minimum variant AC to run callset frequency filters
    ac_min: ClassVar[int] = config_retrieve(['ValidateMOI', 'min_callset_ac_to_filter'])
    ac_threshold: ClassVar[int] = config_retrieve(['ValidateMOI', 'dominant_callset_max_ac'])
    small_af: ClassVar[float] = config_retrieve(['ValidateMOI', 'dominant_callset_max_af'])
    sv_af: ClassVar[float] = config_retrieve(['ValidateMOI', 'dominant_callset_sv_max_af'])

    # a lookup of the attribute name vs. the corresponding configurable filter
    small_dict: ClassVar[dict[str, float | int]] = {
        'gnomad_af': config_retrieve(['ValidateMOI', 'dominant_gnomad_max_af']),
        'gnomad_ac': config_retrieve(['ValidateMOI', 'dominant_gnomad_max_ac']),
        'gnomad_homalt': config_retrieve(['ValidateMOI', 'dominant_gnomad_max_homozygotes']),
    }

    # specific to SVs
    sv_dict: ClassVar[dict[str, float]] = {
        'gnomad_v2.1_sv_AF': config_retrieve(['ValidateMOI', 'dominant_gnomad_sv_max_af']),
    }

    def too_common(self, variant: SmallVariant | StructuralVariant) -> bool:
        """
        Check if a variant is too common in the population

        Args:
            variant (SmallVariant | StructuralVariant): the variant to check

        Returns:
            bool: True if the variant is too common
        """

        # check against each small-variant filter
        if isinstance(variant, SmallVariant):
            for key, threshold in self.small_dict.items():
                if key in variant.info and variant.info[key] > threshold:
                    return True
            if variant.info['ac'] > self.ac_threshold and variant.info['af'] > self.small_af:
                return True

        elif isinstance(variant, StructuralVariant):
            for key, threshold in self.sv_dict.items():
                if key in variant.info and variant.info[key] > threshold:
                    return True
            if variant.info['ac'] > self.ac_min and (
                (variant.info['af'] > self.sv_af) or (variant.info['ac'] > self.ac_threshold)
            ):
                return True

        else:
            raise ValueError('Variant type not recognised')

        return False


@dataclass
class ClinVarFilter:
    """
    This will apply more lenient filters to ClinVar Pathogenic variants
    """

    # minimum variant AC to run callset frequency filters
    ac_threshold: ClassVar[int] = config_retrieve(['ValidateMOI', 'min_callset_ac_to_filter'])

    # a lookup of the attribute name vs. the corresponding configurable filter
    small_dict: ClassVar[dict[str, float]] = {
        'gnomad_af': config_retrieve(['ValidateMOI', 'clinvar_gnomad_max_af']),
    }
    small_af: ClassVar[float] = config_retrieve(['ValidateMOI', 'clinvar_callset_max_af'])

    def too_common(self, variant: SmallVariant) -> bool:
        """
        Check if a variant is too common in the population

        Args:
            variant (SmallVariant): the variant to check

        Returns:
            bool: True if the variant is too common
        """
        for key, threshold in self.small_dict.items():
            if key in variant.info and variant.info[key] > threshold:
                return True

        return variant.info['ac'] > self.ac_threshold and variant.info['af'] > self.small_af


@dataclass
class ClinVarDominantFilter:
    """
    This will apply more lenient filters to ClinVar Pathogenic variants
    Designed to run on Dominant variants
    """

    # minimum variant AC to run callset frequency filters
    ac_threshold: ClassVar[int] = config_retrieve(['ValidateMOI', 'min_callset_ac_to_filter'])

    # a lookup of the attribute name vs. the corresponding configurable filter
    small_dict: ClassVar[dict[str, float]] = {
        'gnomad_af': config_retrieve(['ValidateMOI', 'clinvar_dominant_gnomad_max_af']),
    }
    small_af: ClassVar[float] = config_retrieve(['ValidateMOI', 'clinvar_dominant_callset_max_af'])

    def too_common(self, variant: SmallVariant) -> bool:
        """
        Check if a variant is too common in the population

        Args:
            variant (SmallVariant): the variant to check

        Returns:
            bool: True if the variant is too common
        """
        for key, threshold in self.small_dict.items():
            if key in variant.info and variant.info[key] > threshold:
                return True

        return variant.info['ac'] > self.ac_threshold and variant.info['af'] > self.small_af


class MOIRunner:
    """
    The abstract class for a single MOI runner
    This will be instantiated once per MOI, and run once per related gene, on the collection of all variants in the gene
    """

    def __init__(self, pedigree: PedigreeParser, target_moi: str):
        """
        for each possible MOI, choose the appropriate filters to apply
        ran into a situation where the ID of target_moi didn't match the
        exact same MOI as the IDs were different.

        This logic is only called once per MOI, not once per variant
        """

        self.filter_list: list[BaseMoi] = []

        # for unknown, we catch all possible options?
        # should we be doing both checks for Monoallelic?
        if target_moi == 'Monoallelic':
            self.filter_list = [DominantAutosomal(pedigree=pedigree)]
        elif target_moi in ['Mono_And_Biallelic', 'Unknown']:
            self.filter_list = [
                DominantAutosomal(pedigree=pedigree),
                RecessiveAutosomalHomo(pedigree=pedigree),
                RecessiveAutosomalCH(pedigree=pedigree),
            ]
        elif target_moi == 'Biallelic':
            self.filter_list = [RecessiveAutosomalHomo(pedigree=pedigree), RecessiveAutosomalCH(pedigree=pedigree)]

        elif target_moi == 'Hemi_Mono_In_Female':
            self.filter_list = [XRecessiveMale(pedigree=pedigree), XDominant(pedigree=pedigree)]

        elif target_moi == 'Hemi_Bi_In_Female':
            self.filter_list = [
                XRecessiveMale(pedigree=pedigree),
                XRecessiveFemaleHom(pedigree=pedigree),
                XRecessiveFemaleCH(pedigree=pedigree),
                XPseudoDominantFemale(pedigree=pedigree),
            ]

        else:
            raise KeyError(f'MOI type {target_moi} is not addressed in MOI')

    def run(self, principal_var, comp_het: CompHetDict | None = None, partial_pen: bool = False) -> list[ReportVariant]:
        """
        run method - triggers each relevant inheritance model

        Args:
            principal_var (): the variant we are focused on
            comp_het ():
            partial_pen ():
        """

        if comp_het is None:
            comp_het = {}

        moi_matched = []
        for model in self.filter_list:
            moi_matched.extend(model.run(principal=principal_var, comp_het=comp_het, partial_pen=partial_pen))
        return moi_matched


class BaseMoi:
    """
    Definition of the MOI base class
    """

    def __init__(self, pedigree: PedigreeParser, applied_moi: str):
        """
        base class
        """
        if applied_moi is None:
            raise ValueError('An applied MOI needs to reach the Base Class')
        self.pedigree = pedigree
        self.applied_moi = applied_moi
        self.minimum_alt_depth = config_retrieve(['RunHailFiltering', 'min_alt_depth'], 5)
        self.minimum_depth = config_retrieve(['RunHailFiltering', 'minimum_depth'], 10)
        self.global_filter = GlobalFilter()
        self.clinvar_filter = ClinVarFilter()

    @abstractmethod
    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,
    ) -> list[ReportVariant]:
        """
        run all applicable inheritance patterns and finds good fits
        """

    def single_variant_explains_disease_in_family(
        self,
        sample_id: str,
        called_variants: set[str],
        partial_pen: bool = False,
    ) -> bool:
        """
        check for single variant inheritance
        - find the family ID from this sample
        - iterate through all family members, and check that MOI holds for all

        Return False if a participant fails tests, True if all checked are ok

        NOTE: this called_variants pool is prepared before calling this method.
        If we only want to check hom calls, only send hom calls. If we are checking for
        dominant conditions, we would bundle hom and het calls into this set

        Sex chrom MOI is being broken down in a more granular way, so this may need
        revisiting
        """

        this_member = self.pedigree.participants[sample_id]

        # check for valid inheritance within the immediate trio, if possible
        for member_id in [this_member.father_id, this_member.mother_id]:
            if member_id is None or (iter_member := self.pedigree.participants.get(member_id)) is None:
                continue

            # complete & incomplete penetrance - affected samples must have the variant
            # complete pen. requires participants to be affected if they have the var
            # if any of these combinations occur, fail the family
            if (iter_member.is_affected and member_id not in called_variants) or (
                member_id in called_variants and not partial_pen and iter_member.is_not_affected
            ):
                return False

        return True

    def variant_too_common(self, variant: VARIANT_MODELS) -> bool:
        """
        Check if a variant is too common in the population or callset

        Args:
            variant (VARIANT_MODELS): the variant to check

        Returns:
            bool: True if the variant is too common
        """
        if variant.info.get('categoryboolean1'):
            return self.clinvar_filter.too_common(variant=variant)
        return self.global_filter.too_common(variant)

    def get_family_genotypes(self, variant: VARIANT_MODELS, sample_id: str) -> dict[str, str]:
        """

        Args:
            variant (VARIANT_MODELS):
            sample_id (str): the sample ID to gather genotypes for

        Returns:
            list[str]: a list of all the participants, and GTs
        """

        def get_sample_genotype(member_id: str, sex: int) -> str:
            """
            for this specific member, find the genotype
            Args:
                member_id (str): sample ID in the pedigree
                sex (int): male/female/unknown

            Returns:
                str: text representation of this genotype
            """

            if variant.coordinates.chrom in X_CHROMOSOME:
                if sex == 1 and (member_id in variant.het_samples or member_id in variant.hom_samples):
                    return 'Hemi'

                if member_id in variant.het_samples:
                    return 'Het'
                if member_id in variant.hom_samples:
                    return 'Hom'

            elif member_id in variant.het_samples:
                return 'Het'
            elif member_id in variant.hom_samples:
                return 'Hom'

            return 'WT'

        sample_family_id = self.pedigree.participants[sample_id].family_id
        return {
            member.sample_id: get_sample_genotype(member_id=member.sample_id, sex=member.sex)
            for member in self.pedigree.by_family[sample_family_id]
        }

    def comp_het_explains_disease_in_family(
        self,
        sample_id: str,
        variant_1: VARIANT_MODELS,
        variant_2: VARIANT_MODELS,
    ) -> bool:
        """
        use parents to accept or dismiss the comp-het
        If the 'comp-het' pair are inherited from a single parent, they are in cis
        rather than trans, and reporting as a comp-het would be misleading.

        compound het is inherently not inherited from a single parent, so rule out
        when either parent has both, or either parent is affected

        Args:
            sample_id (str): sample ID to check for
            variant_1 (VARIANT_MODELS): first variant of comp-het pair
            variant_2 (VARIANT_MODELS): second variant of comp-het pair

        Returns:
            bool: True if these two variants form a comp-het
        """

        # if both vars are present in a single parent: not a compound het
        # or if the parent is affected: not causative
        sample_ped_entry = self.pedigree.participants[sample_id]
        for parent in [sample_ped_entry.mother_id, sample_ped_entry.father_id]:
            # skip to prevent crashing on !trios, or participants without a pedigree row of their own.
            if parent is None or self.pedigree.participants.get(parent, None) is None:
                continue

            if ((parent in variant_1.het_samples) and (parent in variant_2.het_samples)) or self.pedigree.participants[
                parent
            ].affected == 2:
                return False
        return True


class DominantAutosomal(BaseMoi):
    """This class can also be called by the X-linked Dominant, in which case the Applied_MOI by name is overridden."""

    def __init__(self, pedigree: PedigreeParser, applied_moi: str = 'Autosomal Dominant'):
        """
        Simplest: AD MOI
        """

        super().__init__(pedigree=pedigree, applied_moi=applied_moi)
        self.global_filter = DominantFilter()
        self.clinvar_filter = ClinVarDominantFilter()

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,  # noqa: ARG002
        partial_pen: bool = False,
    ) -> list[ReportVariant]:
        """Simplest MOI, exclusions based on HOM count and AF."""

        classifications = []

        if self.variant_too_common(principal):
            return classifications

        # autosomal dominant doesn't require support, but consider het and hom
        samples_with_this_variant = principal.het_samples.union(principal.hom_samples)
        for sample_id in samples_with_this_variant:
            # skip primary analysis for unaffected members
            # we require this specific sample to be categorised
            # force a minimum depth on the proband call
            # and a minimum number of alt reads supporting
            if (
                self.pedigree.participants[sample_id].is_not_affected
                or not principal.sample_category_check(sample_id, allow_support=False)
                or principal.insufficient_read_depth(
                    sample_id,
                    self.minimum_depth,
                    var_is_cat_1=principal.info.get('categoryboolean1'),
                )
                or principal.insufficient_alt_depth(sample_id, self.minimum_alt_depth)
            ):
                continue

            # check if this is a candidate for dominant inheritance
            if self.single_variant_explains_disease_in_family(
                sample_id=sample_id,
                called_variants=samples_with_this_variant,
                partial_pen=partial_pen,
            ):
                classifications.append(
                    ReportVariant(
                        sample=sample_id,
                        family=self.pedigree.participants[sample_id].family_id,
                        gene=principal.info.get('gene_id'),
                        var_data=principal,
                        categories=principal.category_values(sample_id),
                        reasons=self.applied_moi,
                        genotypes=self.get_family_genotypes(variant=principal, sample_id=sample_id),
                        flags=principal.get_sample_flags(sample_id),
                        independent=True,
                    ),
                )

        return classifications


class RecessiveAutosomalCH(BaseMoi):
    """
    inheritance test class for Recessive inheritance
    requires single hom variant, or compound het
    """

    def __init__(self, pedigree: PedigreeParser, applied_moi: str = 'Autosomal Recessive Comp-Het'):
        """ """
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,  # noqa: ARG002
    ) -> list[ReportVariant]:
        """
        Valid if present as compound het
        counts as being phased if a compound het is split between parents
        """

        if comp_het is None:
            comp_het = {}

        classifications = []

        # reject support for dominant MOI, apply checks based on var type
        if self.variant_too_common(principal):
            return classifications

        # if hets are present, try and find support
        for sample_id in principal.het_samples:
            # skip primary analysis for unaffected members
            # this sample must be categorised - check Cat 4 contents
            if (
                self.pedigree.participants[sample_id].affected != 2
                or (not principal.sample_category_check(sample_id, allow_support=True))
                or principal.insufficient_read_depth(
                    sample_id,
                    self.minimum_depth,
                    principal.info.get('categoryboolean1'),
                )
                or principal.insufficient_alt_depth(sample_id, self.minimum_alt_depth)
            ):
                continue

            for partner in comp_het[sample_id].get(principal.coordinates.string_format, []):
                if (
                    partner.insufficient_read_depth(
                        sample_id,
                        self.minimum_depth,
                        partner.info.get('categoryboolean1'),
                    )
                    or self.variant_too_common(partner)
                    or partner.insufficient_alt_depth(sample_id, self.minimum_alt_depth)
                    or not any(
                        [
                            principal.sample_category_check(sample_id, allow_support=False),
                            partner.sample_category_check(sample_id, allow_support=False),
                        ],
                    )
                ):
                    continue

                # check if this is a candidate for comp-het inheritance
                if self.comp_het_explains_disease_in_family(
                    sample_id=sample_id,
                    variant_1=principal,
                    variant_2=partner,
                ):
                    classifications.append(
                        ReportVariant(
                            sample=sample_id,
                            family=self.pedigree.participants[sample_id].family_id,
                            gene=principal.info.get('gene_id'),
                            var_data=principal,
                            categories=principal.category_values(sample_id),
                            reasons=self.applied_moi,
                            genotypes=self.get_family_genotypes(variant=principal, sample_id=sample_id),
                            support_vars={partner.info['var_link']},
                            flags=principal.get_sample_flags(sample_id) | partner.get_sample_flags(sample_id),
                            independent=False,
                        ),
                    )

        return classifications


class RecessiveAutosomalHomo(BaseMoi):
    """
    inheritance test class for Recessive inheritance
    requires single hom variant, or compound het
    """

    def __init__(self, pedigree: PedigreeParser, applied_moi: str = 'Autosomal Recessive Homozygous'):
        """ """
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,  # noqa: ARG002
        partial_pen: bool = False,
    ) -> list[ReportVariant]:
        """
        explicitly tests HOMs

        Args:
            principal (VARIANT_MODELS): main variant being evaluated
            comp_het (dict): comp-het partners
            partial_pen (bool):

        Returns:
            list[ReportVariant]: data object if RecessiveAutosomal fits
        """

        classifications = []

        # remove if too many homs are present in population databases
        if self.variant_too_common(principal):
            return classifications

        for sample_id in principal.hom_samples:
            # skip primary analysis for unaffected members
            # require this sample to be categorised - check Sample contents
            # minimum depth of call
            if (
                self.pedigree.participants[sample_id].affected != 2
                or not (principal.sample_category_check(sample_id, allow_support=False))
                or principal.insufficient_read_depth(
                    sample_id,
                    self.minimum_depth,
                    principal.info.get('categoryboolean1'),
                )
                or principal.insufficient_alt_depth(sample_id, self.minimum_alt_depth)
            ):
                continue

            # check if this is a possible candidate for homozygous inheritance
            if self.single_variant_explains_disease_in_family(
                sample_id=sample_id,
                called_variants=principal.hom_samples,
                partial_pen=partial_pen,
            ):
                classifications.append(
                    ReportVariant(
                        sample=sample_id,
                        family=self.pedigree.participants[sample_id].family_id,
                        gene=principal.info.get('gene_id'),
                        var_data=principal,
                        categories=principal.category_values(sample_id),
                        genotypes=self.get_family_genotypes(variant=principal, sample_id=sample_id),
                        reasons=self.applied_moi,
                        flags=principal.get_sample_flags(sample_id),
                        independent=True,
                    ),
                )

        return classifications


class XDominant(BaseMoi):
    """
    for males and females, accept het
    effectively the same as DominantAutosomal?
    just override the type, and use AD

    GATK MALES ARE CALLED HOM (OR HET pseudo-autosomal?)
    re-implement here, but don't permit Male X-Homs
    """

    def __init__(self, pedigree: PedigreeParser, applied_moi: str = 'X_Dominant'):
        """accept male hets and homs, and female hets without support."""
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)
        self.global_filter = DominantFilter()
        self.clinvar_filter = ClinVarDominantFilter()

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,  # noqa: ARG002
        partial_pen: bool = False,
    ) -> list[ReportVariant]:
        """
        If variant is present and sufficiently rare, we take it. discarded if support only.
        """

        classifications = []

        # more stringent Pop.Freq checks for dominant - hemi restriction
        if self.variant_too_common(principal):
            return classifications

        # all samples which have a variant call
        samples_with_this_variant = principal.het_samples.union(principal.hom_samples)

        for sample_id in samples_with_this_variant:
            # skip primary analysis for unaffected members
            # we require this specific sample to be categorised (non-support)
            # force minimum depth
            if (
                self.pedigree.participants[sample_id].affected != 2
                or not (principal.sample_category_check(sample_id, allow_support=False))
                or principal.insufficient_read_depth(
                    sample_id,
                    self.minimum_depth,
                    principal.info.get('categoryboolean1'),
                )
                or principal.insufficient_alt_depth(sample_id, self.minimum_alt_depth)
            ):
                continue

            # check if this is a candidate for dominant inheritance
            if self.single_variant_explains_disease_in_family(
                sample_id=sample_id,
                called_variants=samples_with_this_variant,
                partial_pen=partial_pen,
            ):
                classifications.append(
                    ReportVariant(
                        sample=sample_id,
                        family=self.pedigree.participants[sample_id].family_id,
                        gene=principal.info.get('gene_id'),
                        var_data=principal,
                        categories=principal.category_values(sample_id),
                        reasons=self.applied_moi,
                        genotypes=self.get_family_genotypes(variant=principal, sample_id=sample_id),
                        flags=principal.get_sample_flags(sample_id),
                        independent=True,
                    ),
                )
        return classifications


class XPseudoDominantFemale(BaseMoi):
    """
    X-Dominant method which only evaluates females, on the basis that a healthy allele could be inactivated
    A special case of X-Dominant, where we consider Het. variants, used to assess genes which are on X and in genes
    associated with a Biallelic MOI

    Basically a Dominant MOI to be applied to Recessive genes, and results will be labelled as cautionary
    """

    def __init__(self, pedigree: PedigreeParser, applied_moi: str = 'X_PseudoDominant'):
        """Accept male hets and homs, and female hets without support/"""
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)
        self.global_filter = DominantFilter()
        self.clinvar_filter = ClinVarDominantFilter()

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,  # noqa: ARG002
        partial_pen: bool = False,
    ) -> list[ReportVariant]:
        """If variant is present and sufficiently rare, we take it. Discarded if support."""

        # unused in this class, we always run this with partial penetrance
        _unused = partial_pen

        classifications = []

        # more stringent Pop.Freq checks for dominant - hemi restriction
        if self.variant_too_common(principal):
            return classifications

        # all females which have a variant call
        considered_females = {sam for sam in principal.het_samples if self.pedigree.participants[sam].is_female}
        all_with_variant = principal.het_samples.union(principal.hom_samples)
        for sample_id in considered_females:
            # skip primary analysis for unaffected members
            # we require this specific sample to be categorised
            # force minimum depth
            if (
                self.pedigree.participants[sample_id].is_not_affected
                or not (principal.sample_category_check(sample_id, allow_support=False))
                or principal.insufficient_read_depth(
                    sample_id,
                    self.minimum_depth,
                    principal.info.get('categoryboolean1'),
                )
                or principal.insufficient_alt_depth(sample_id, self.minimum_alt_depth)
            ):
                continue

            # check if this is a candidate for dominant inheritance
            # as we're allowing for flexible 'penetrance' in females, we send all het and hom variants, but allow for a
            # partial penetrance check - the participants can have the variant, but not be affected. They may not be
            # affected without the variant call.
            # There's a slight breakdown here as the males should be interpreted under a full penetrance model, and
            # females under partial penetrance, but that's not trivial without creating a second familial check method.
            # Leaving that aside now as the current implementation is pretty central to the algorithm. Will revisit if
            # this is noisy.
            if self.single_variant_explains_disease_in_family(
                sample_id=sample_id,
                called_variants=all_with_variant,
                partial_pen=True,
            ):
                classifications.append(
                    ReportVariant(
                        sample=sample_id,
                        family=self.pedigree.participants[sample_id].family_id,
                        gene=principal.info.get('gene_id'),
                        var_data=principal,
                        categories=principal.category_values(sample_id),
                        reasons=self.applied_moi,
                        genotypes=self.get_family_genotypes(variant=principal, sample_id=sample_id),
                        flags=principal.get_sample_flags(sample_id)
                        | {'Affected female with heterozygous variant in XLR gene'},
                        independent=True,
                    ),
                )
        return classifications


class XRecessiveMale(BaseMoi):
    """
    accept all male non-ref GTs - male variants HOM because GATK
    effectively the same as AutosomalDominant?
    """

    def __init__(self, pedigree: PedigreeParser, applied_moi: str = 'X_Male'):
        """Set parameters specific to male X tests."""
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)
        self.global_filter = DominantFilter()
        self.clinvar_filter = ClinVarDominantFilter()

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,  # noqa: ARG002
        partial_pen: bool = False,
    ) -> list[ReportVariant]:
        """Apply filtering specific to variants on X in Males."""

        classifications = []

        # remove from analysis if too many homs are present in population databases
        if self.variant_too_common(principal):
            return classifications

        # combine het and hom here, we don't trust the variant callers
        # if hemi count is too high, don't consider males
        # never consider support variants on X for males
        males = {
            sam for sam in principal.het_samples.union(principal.hom_samples) if self.pedigree.participants[sam].is_male
        }

        for sample_id in males:
            # specific affected sample category check, never consider support on X for males
            if (
                self.pedigree.participants[sample_id].is_not_affected
                or not (principal.sample_category_check(sample_id, allow_support=False))
                or principal.insufficient_read_depth(
                    sample_id,
                    self.minimum_depth,
                    principal.info.get('categoryboolean1'),
                )
                or principal.insufficient_alt_depth(sample_id, self.minimum_alt_depth)
            ):
                continue

            # check if this is a possible candidate for homozygous inheritance
            if self.single_variant_explains_disease_in_family(
                sample_id=sample_id,
                called_variants=males,
                partial_pen=partial_pen,
            ):
                classifications.append(
                    ReportVariant(
                        sample=sample_id,
                        family=self.pedigree.participants[sample_id].family_id,
                        gene=principal.info.get('gene_id'),
                        var_data=principal,
                        categories=principal.category_values(sample_id),
                        genotypes=self.get_family_genotypes(variant=principal, sample_id=sample_id),
                        reasons=self.applied_moi,
                        flags=principal.get_sample_flags(sample_id),
                        independent=True,
                    ),
                )
        return classifications


class XRecessiveFemaleHom(BaseMoi):
    """
    only consider HOM females
    """

    def __init__(self, pedigree: PedigreeParser, applied_moi: str = 'X_Recessive HOM Female'):
        """Set parameters specific to recessive tests."""
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,  # noqa: ARG002
        partial_pen: bool = False,
    ) -> list[ReportVariant]:
        """ """

        classifications = []

        # remove from analysis if too many homs are present in population databases
        if self.variant_too_common(principal):
            return classifications

        # never consider support homs
        samples_to_check = {sam for sam in principal.hom_samples if self.pedigree.participants[sam].is_female}

        for sample_id in samples_to_check:
            # specific affected sample category check
            if (
                self.pedigree.participants[sample_id].is_not_affected
                or not principal.sample_category_check(sample_id, allow_support=False)
                or principal.insufficient_read_depth(
                    sample_id,
                    self.minimum_depth,
                    principal.info.get('categoryboolean1'),
                )
                or principal.insufficient_alt_depth(sample_id, self.minimum_alt_depth)
            ):
                continue

            # check if this is a possible candidate for homozygous inheritance
            if self.single_variant_explains_disease_in_family(
                sample_id=sample_id,
                called_variants=samples_to_check,
                partial_pen=partial_pen,
            ):
                classifications.append(
                    ReportVariant(
                        sample=sample_id,
                        family=self.pedigree.participants[sample_id].family_id,
                        gene=principal.info.get('gene_id'),
                        var_data=principal,
                        categories=principal.category_values(sample_id),
                        genotypes=self.get_family_genotypes(variant=principal, sample_id=sample_id),
                        reasons=self.applied_moi,
                        flags=principal.get_sample_flags(sample_id),
                        independent=True,
                    ),
                )
        return classifications


class XRecessiveFemaleCH(BaseMoi):
    """
    ignore males, accept female comp-het only
    """

    def __init__(self, pedigree: PedigreeParser, applied_moi: str = 'X_RecessiveFemaleCompHet'):
        """Set parameters specific to recessive tests on X in females, specific to Compound-Heterozygous pairs."""
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,  # noqa: ARG002
    ) -> list[ReportVariant]:
        """ """

        if comp_het is None:
            comp_het = {}

        classifications = []

        # remove from analysis if too many homs are present in population databases
        if self.variant_too_common(principal):
            return classifications
        het_females = {sam for sam in principal.het_samples if self.pedigree.participants[sam].is_female}

        # if het females are present, try and find support
        for sample_id in het_females:
            # don't run primary analysis for unaffected
            # we require this specific sample to be categorised - check Cat 4 contents
            if (
                self.pedigree.participants[sample_id].is_not_affected
                or not principal.sample_category_check(sample_id, allow_support=True)
                or principal.insufficient_read_depth(
                    sample_id,
                    self.minimum_depth,
                    principal.info.get('categoryboolean1'),
                )
                or principal.insufficient_alt_depth(sample_id, self.minimum_alt_depth)
            ):
                continue

            for partner in comp_het[sample_id].get(principal.coordinates.string_format, []):
                # allow for de novo check - also screen out high-AF partners
                # check for minimum depth and alt support in partner
                if (
                    self.variant_too_common(partner)
                    or partner.insufficient_alt_depth(sample_id, self.minimum_alt_depth)
                    or partner.insufficient_read_depth(
                        sample_id,
                        self.minimum_depth,
                        partner.info.get('categoryboolean1'),
                    )
                    or not any(
                        [
                            principal.sample_category_check(sample_id, allow_support=False),
                            partner.sample_category_check(sample_id, allow_support=False),
                        ],
                    )
                ):
                    continue

                if self.comp_het_explains_disease_in_family(
                    sample_id=sample_id,
                    variant_1=principal,
                    variant_2=partner,
                ):
                    classifications.append(
                        ReportVariant(
                            sample=sample_id,
                            family=self.pedigree.participants[sample_id].family_id,
                            gene=principal.info.get('gene_id'),
                            var_data=principal,
                            categories=principal.category_values(sample_id),
                            reasons=self.applied_moi,
                            genotypes=self.get_family_genotypes(variant=principal, sample_id=sample_id),
                            # needs to comply with Seqr
                            support_vars={partner.info['var_link']},
                            flags=principal.get_sample_flags(sample_id) | partner.get_sample_flags(sample_id),
                            independent=False,
                        ),
                    )

        return classifications
