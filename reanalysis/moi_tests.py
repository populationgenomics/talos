"""
A number of classes, each representing one Mode of Inheritance
One class (MoiRunner) to run all the appropriate MOIs on a variant

Reduce the PanelApp plain text MOI description into a few categories
We then run a permissive MOI match for the variant
"""
# mypy: ignore-errors
from abc import abstractmethod

from peddy.peddy import Ped, PHENOTYPE

from cpg_utils.config import get_config

from reanalysis.utils import (
    AbstractVariant,
    CompHetDict,
    MinimalVariant,
    ReportedVariant,
    X_CHROMOSOME,
)

# config keys to use for dominant MOI tests
GNOMAD_RARE_THRESHOLD = 'gnomad_dominant'
GNOMAD_AD_AC_THRESHOLD = 'gnomad_max_ac_dominant'
GNOMAD_DOM_HOM_THRESHOLD = 'gnomad_max_homs_dominant'
GNOMAD_REC_HOM_THRESHOLD = 'gnomad_max_homs_recessive'
GNOMAD_HEMI_THRESHOLD = 'gnomad_max_hemi'
INFO_HOMS = {'gnomad_hom', 'gnomad_ex_hom'}
INFO_HEMI = {'gnomad_hemi', 'gnomad_ex_hemi'}
PEDDY_AFFECTED = PHENOTYPE().AFFECTED


def check_for_second_hit(
    first_variant: str, comp_hets: CompHetDict, sample: str
) -> list[AbstractVariant]:
    """
    checks for a second hit partner in this gene

    Example formatting of the comp-het dict
    {
        "SampleID": {
            "12-52287177-T-C": [
                AbstractVariant(12-52287180-TGG-T)
            ],
            "12-52287180-TGG-T": [
                AbstractVariant(12-52287177-T-C)
            ]
        } ...
    }

    :param first_variant: string representation of variant1
    :param comp_hets: lookup dict for compound hets
    :param sample: ID string
    :return:
    """

    # check if the sample has any comp-hets
    if sample not in comp_hets.keys():
        return []

    sample_dict = comp_hets.get(sample)
    return sample_dict.get(first_variant, [])


class MOIRunner:
    """
    pass
    """

    def __init__(self, pedigree: Ped, target_moi: str):
        """
        for each possible MOI, choose the appropriate filters to apply
        ran into a situation where the ID of target_moi didn't match the
        exact same MOI as the IDs were different.

        This logic is only called once per MOI, not once per variant

        :param pedigree:
        :param target_moi:
        """

        # for unknown, we catch all possible options?
        # should we be doing both checks for Monoallelic?
        if target_moi == 'Monoallelic':
            self.filter_list = [
                DominantAutosomal(pedigree=pedigree),
            ]
        elif target_moi in ['Mono_And_Biallelic', 'Unknown']:
            self.filter_list = [
                DominantAutosomal(pedigree=pedigree),
                RecessiveAutosomalHomo(pedigree=pedigree),
                RecessiveAutosomalCH(pedigree=pedigree),
            ]
        elif target_moi == 'Biallelic':
            self.filter_list = [
                RecessiveAutosomalHomo(pedigree=pedigree),
                RecessiveAutosomalCH(pedigree=pedigree),
            ]

        elif target_moi == 'Hemi_Mono_In_Female':
            self.filter_list = [
                XRecessiveMale(pedigree=pedigree),
                XDominant(pedigree=pedigree),
            ]

        elif target_moi == 'Hemi_Bi_In_Female':
            self.filter_list = [
                XRecessiveMale(pedigree=pedigree),
                XRecessiveFemaleHom(pedigree=pedigree),
                XRecessiveFemaleCH(pedigree=pedigree),
            ]

        else:
            raise KeyError(f'MOI type {target_moi} is not addressed in MOI')

    def run(
        self,
        principal_var,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,
    ) -> list[ReportedVariant]:
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
            moi_matched.extend(
                model.run(
                    principal=principal_var, comp_het=comp_het, partial_pen=partial_pen
                )
            )
        return moi_matched


class BaseMoi:
    """
    Definition of the MOI base class
    """

    def __init__(self, pedigree: Ped, applied_moi: str):
        """
        base class
        """
        if applied_moi is None:
            raise ValueError('An applied MOI needs to reach the Base Class')
        self.pedigree = pedigree
        self.applied_moi = applied_moi
        self.minimum_depth = get_config()['filter'].get('minimum_depth', 10)

    @abstractmethod
    def run(
        self,
        principal: AbstractVariant,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,
    ) -> list[ReportedVariant]:
        """
        run all applicable inheritance patterns and finds good fits
        """

    def check_familial_inheritance(
        self, sample_id: str, called_variants: set[str], partial_pen: bool = False
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

        :param sample_id:
        :param called_variants: the set of sample_ids which have this variant
        :param partial_pen: if True, permit unaffected has variant call
        """

        # iterate through all family members, no interested in directionality
        # of relationships at the moment
        for member in self.pedigree.families[self.pedigree[sample_id].family_id]:

            # complete & incomplete penetrance - affected samples must have the variant
            # complete pen. requires participants to be affected if they have the var
            # if any of these combinations occur, fail the family
            if (
                member.affected == PEDDY_AFFECTED
                and member.sample_id not in called_variants
            ) or (
                member.sample_id in called_variants
                and not partial_pen
                and not member.affected == PEDDY_AFFECTED
            ):
                return False

        return True

    def get_family_genotypes(
        self, variant: AbstractVariant, sample_id: str
    ) -> dict[str, str]:
        """

        Args:
            variant (AbstractVariant):
            sample_id (str): the sample ID to gather genotypes for

        Returns:
            list[str]: a list of all the participants, and GTs
        """

        def get_sample_genotype(member_id: str, sex: str) -> str:
            """
            for this specific member, find the genotype
            Args:
                member_id (str): sample ID in the pedigree
                sex (str): male/female/unknown

            Returns:
                str: text representation of this genotype
            """

            if variant.coords.chrom in X_CHROMOSOME:
                if sex == 'male' and (
                    member_id in variant.het_samples or member_id in variant.hom_samples
                ):
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

        sample_ped_entry = self.pedigree[sample_id]
        family = self.pedigree.families[sample_ped_entry.family_id]
        return {
            member.sample_id: get_sample_genotype(
                member_id=member.sample_id, sex=member.sex
            )
            for member in family.samples
        }

    @staticmethod
    def check_frequency(info: dict, keys: set[str], threshold: int | float) -> bool:
        """
        Method to check multiple info keys against a single threshold
        This just reduces the line count, as this is called a bunch of times

        Args:
            info (): the dict of values for this dict
            keys (): the iterable of keys to check
            threshold (): the threshold to test against

        Returns:
            True if any of the info attributes is above the threshold
        """
        return any({info.get(key, 0) > threshold for key in keys})

    def check_comp_het(
        self, sample_id: str, variant_1: AbstractVariant, variant_2: AbstractVariant
    ) -> bool:
        """
        use parents to accept or dismiss the comp-het
        If the 'comp-het' pair are inherited from a single parent, they are in cis
        rather than trans, and reporting as a comp-het would be misleading.

        compound het is inherently not inherited from a single parent, so rule out
        when either parent has both, or either parent is affected

        Args:
            sample_id (str): sample ID to check for
            variant_1 (AbstractVariant): first variant of comp-het pair
            variant_2 (AbstractVariant): second variant of comp-het pair

        Returns:
            bool: True if these two variants form a comp-het
        """

        # if both vars are present in a single parent: not a compound het
        # or if the parent is affected: not causative
        sample_ped_entry = self.pedigree[sample_id]
        for parent in [sample_ped_entry.mom, sample_ped_entry.dad]:
            # skip to prevent crashing on !trios
            if parent is None:
                continue

            if (
                (parent.sample_id in variant_1.het_samples)
                and (parent.sample_id in variant_2.het_samples)
            ) or parent.affected == PEDDY_AFFECTED:
                return False
        return True


class DominantAutosomal(BaseMoi):
    """
    This class can also be called by the X-linked Dominant, in which case the
    Applied_MOI by name is overridden
    """

    def __init__(
        self,
        pedigree: Ped,
        applied_moi: str = 'Autosomal Dominant',
    ):
        """
        Simplest: AD MOI
        """
        self.ad_threshold = get_config()['moi_tests'][GNOMAD_RARE_THRESHOLD]
        self.ac_threshold = get_config()['moi_tests'][GNOMAD_AD_AC_THRESHOLD]
        self.hom_threshold = get_config()['moi_tests'][GNOMAD_DOM_HOM_THRESHOLD]
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: AbstractVariant,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,
    ) -> list[ReportedVariant]:
        """
        Simplest MOI, exclusions based on HOM count and AF
        Args:
            principal ():
            comp_het ():
            partial_pen ():
        """

        classifications = []

        # more stringent Pop.Freq checks for dominant
        # reject support for dominant MOI
        if principal.support_only or (
            (
                self.check_frequency(principal.info, INFO_HOMS, self.hom_threshold)
                or self.check_frequency(
                    principal.info, {'gnomad_ac'}, self.ac_threshold
                )
                or self.check_frequency(
                    principal.info, {'gnomad_af'}, self.ad_threshold
                )
            )
        ):
            return classifications

        # autosomal dominant doesn't require support, but consider het and hom
        samples_with_this_variant = principal.het_samples.union(principal.hom_samples)

        for sample_id in samples_with_this_variant:

            # skip primary analysis for unaffected members
            # we require this specific sample to be categorised
            # force a minimum depth on the proband call
            if not (
                self.pedigree[sample_id].affected == PEDDY_AFFECTED
                and principal.sample_category_check(sample_id, allow_support=False)
            ) or (
                (
                    principal.depths[sample_id] < self.minimum_depth
                    and not principal.info.get('categoryboolean1')
                )
                and 'SVTYPE' not in principal.info
            ):
                continue

            # check if this is a candidate for dominant inheritance
            if not self.check_familial_inheritance(
                sample_id=sample_id,
                called_variants=samples_with_this_variant,
                partial_pen=partial_pen,
            ):
                continue

            classifications.append(
                ReportedVariant(
                    sample=sample_id,
                    family=self.pedigree[sample_id].family_id,
                    gene=principal.info.get('gene_id'),
                    var_data=MinimalVariant(variant=principal, sample=sample_id),
                    reasons={self.applied_moi},
                    genotypes=self.get_family_genotypes(
                        variant=principal, sample_id=sample_id
                    ),
                    flags=principal.get_sample_flags(sample_id),
                )
            )

        return classifications


class RecessiveAutosomalCH(BaseMoi):
    """
    inheritance test class for Recessive inheritance
    requires single hom variant, or compound het
    """

    def __init__(
        self,
        pedigree: Ped,
        applied_moi: str = 'Autosomal Recessive Comp-Het',
    ):
        """ """
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: AbstractVariant,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,
    ) -> list[ReportedVariant]:
        """
        valid if present as compound het
        counts as being phased if a compound het is split between parents

        Args:
            principal (AbstractVariant): main variant being evaluated
            comp_het (dict): comp-het partners
            partial_pen (bool):

        Returns:
            list[ReportedVariant]: data object if RecessiveAutosomal fits
        """

        if comp_het is None:
            comp_het = {}

        classifications = []

        # if hets are present, try and find support
        for sample_id in principal.het_samples:

            # skip primary analysis for unaffected members
            # this sample must be categorised - check Cat 4 contents
            if (
                not (
                    self.pedigree[sample_id].affected == PEDDY_AFFECTED
                    and principal.sample_category_check(sample_id, allow_support=True)
                )
            ) or (
                principal.depths[sample_id] < self.minimum_depth
                and not principal.info.get('categoryboolean1')
            ):
                continue

            for partner_variant in check_for_second_hit(
                first_variant=principal.coords.string_format,
                comp_hets=comp_het,
                sample=sample_id,
            ):

                # skip the double-support scenario
                if (
                    principal.sample_support_only(sample_id)
                    and partner_variant.sample_support_only(sample_id)
                ) or (
                    partner_variant.depths[sample_id] < self.minimum_depth
                    and not partner_variant.info.get('categoryboolean1')
                ):
                    continue

                # categorised for this specific sample, allow support in partner
                # - also screen out high-AF partners
                if not partner_variant.sample_category_check(
                    sample_id, allow_support=True
                ):
                    continue

                # check if this is a candidate for comp-het inheritance
                if not self.check_comp_het(
                    sample_id=sample_id, variant_1=principal, variant_2=partner_variant
                ):
                    continue
                classifications.append(
                    ReportedVariant(
                        sample=sample_id,
                        family=self.pedigree[sample_id].family_id,
                        gene=principal.info.get('gene_id'),
                        var_data=MinimalVariant(principal, sample_id),
                        reasons={self.applied_moi},
                        genotypes=self.get_family_genotypes(
                            variant=principal, sample_id=sample_id
                        ),
                        support_vars={partner_variant.coords.string_format},
                        flags=principal.get_sample_flags(sample_id)
                        + partner_variant.get_sample_flags(sample_id),
                    ),
                )

        return classifications


class RecessiveAutosomalHomo(BaseMoi):
    """
    inheritance test class for Recessive inheritance
    requires single hom variant, or compound het
    """

    def __init__(
        self,
        pedigree: Ped,
        applied_moi: str = 'Autosomal Recessive Homozygous',
    ):
        """ """
        self.hom_threshold = get_config()['moi_tests'][GNOMAD_REC_HOM_THRESHOLD]
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: AbstractVariant,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,
    ) -> list[ReportedVariant]:
        """
        explicitly tests HOMs

        Args:
            principal (AbstractVariant): main variant being evaluated
            comp_het (dict): comp-het partners
            partial_pen (bool):

        Returns:
            list[ReportedVariant]: data object if RecessiveAutosomal fits
        """

        classifications = []

        # remove if too many homs are present in population databases
        if principal.support_only or (
            self.check_frequency(principal.info, INFO_HOMS, self.hom_threshold)
            and not principal.info.get('categoryboolean1')
        ):
            return classifications

        for sample_id in principal.hom_samples:

            # skip primary analysis for unaffected members
            # require this sample to be categorised - check Sample contents
            # minimum depth of call
            if (
                not (
                    self.pedigree[sample_id].affected == PEDDY_AFFECTED
                    and principal.sample_category_check(sample_id, allow_support=False)
                )
            ) or (
                principal.depths[sample_id] < self.minimum_depth
                and not principal.info.get('categoryboolean1')
            ):
                continue

            # check if this is a possible candidate for homozygous inheritance
            if not self.check_familial_inheritance(
                sample_id=sample_id,
                called_variants=principal.hom_samples,
                partial_pen=partial_pen,
            ):
                continue

            classifications.append(
                ReportedVariant(
                    sample=sample_id,
                    family=self.pedigree[sample_id].family_id,
                    gene=principal.info.get('gene_id'),
                    var_data=MinimalVariant(principal, sample_id),
                    genotypes=self.get_family_genotypes(
                        variant=principal, sample_id=sample_id
                    ),
                    reasons={self.applied_moi},
                    flags=principal.get_sample_flags(sample_id),
                )
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

    def __init__(self, pedigree: Ped, applied_moi: str = 'X_Dominant'):
        """
        accept male hets and homs, and female hets without support
        :param pedigree:
        :param applied_moi:
        """
        self.ad_threshold = get_config()['moi_tests'][GNOMAD_RARE_THRESHOLD]
        self.ac_threshold = get_config()['moi_tests'][GNOMAD_AD_AC_THRESHOLD]
        self.hom_threshold = get_config()['moi_tests'][GNOMAD_DOM_HOM_THRESHOLD]
        self.hemi_threshold = get_config()['moi_tests'][GNOMAD_HEMI_THRESHOLD]
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: AbstractVariant,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,
    ) -> list[ReportedVariant]:
        """
        if variant is present and sufficiently rare, we take it
        discarded if support

        Args:
            principal ():
            comp_het ():
            partial_pen ():
        """

        classifications = []

        if principal.support_only:
            return classifications

        # never apply dominant MOI to support variants
        # more stringent Pop.Freq checks for dominant - hemi restriction
        if (
            self.check_frequency(principal.info, INFO_HOMS, self.hom_threshold)
            or self.check_frequency(principal.info, {'gnomad_ad'}, self.ad_threshold)
            or self.check_frequency(principal.info, {'gnomad_ac'}, self.ac_threshold)
            or self.check_frequency(principal.info, INFO_HEMI, self.hemi_threshold)
        ):
            return classifications

        # all samples which have a variant call
        samples_with_this_variant = principal.het_samples.union(principal.hom_samples)

        for sample_id in samples_with_this_variant:

            # skip primary analysis for unaffected members
            # we require this specific sample to be categorised
            # force minimum depth
            if (
                not (
                    principal.sample_category_check(sample_id, allow_support=False)
                    and self.pedigree[sample_id].affected == PEDDY_AFFECTED
                )
            ) or (
                principal.depths[sample_id] < self.minimum_depth
                and not principal.info.get('categoryboolean1')
            ):
                continue

            # check if this is a candidate for dominant inheritance
            if not self.check_familial_inheritance(
                sample_id=sample_id,
                called_variants=samples_with_this_variant,
                partial_pen=partial_pen,
            ):
                continue

            classifications.append(
                ReportedVariant(
                    sample=sample_id,
                    family=self.pedigree[sample_id].family_id,
                    gene=principal.info.get('gene_id'),
                    var_data=MinimalVariant(principal, sample_id),
                    reasons={self.applied_moi},
                    genotypes=self.get_family_genotypes(
                        variant=principal, sample_id=sample_id
                    ),
                    flags=principal.get_sample_flags(sample_id),
                )
            )
        return classifications


class XRecessiveMale(BaseMoi):
    """
    accept all male non-ref GTs - male variants HOM because GATK
    effectively the same as AutosomalDominant?
    """

    def __init__(
        self,
        pedigree: Ped,
        applied_moi: str = 'X_Male',
    ):
        """
        set parameters specific to male X tests

        Args:
            pedigree ():
            applied_moi ():
        """

        self.hom_dom_threshold = get_config()['moi_tests'][GNOMAD_DOM_HOM_THRESHOLD]
        self.hemi_threshold = get_config()['moi_tests'][GNOMAD_HEMI_THRESHOLD]
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: AbstractVariant,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,
    ) -> list[ReportedVariant]:
        """
        Args:
            principal ():
            comp_het ():
            partial_pen ():
        """

        classifications = []

        # never consider support on X for males
        # remove from analysis if too many homs are present in population databases
        if principal.support_only or (
            (
                self.check_frequency(principal.info, INFO_HOMS, self.hom_dom_threshold)
                or self.check_frequency(principal.info, INFO_HEMI, self.hemi_threshold)
            )
        ):
            return classifications

        # combine het and hom here, we don't trust the variant callers
        # if hemi count is too high, don't consider males
        # never consider support variants on X for males
        males = {
            sam
            for sam in principal.het_samples.union(principal.hom_samples)
            if self.pedigree[sam].sex == 'male'
        }

        for sample_id in males:

            # specific affected sample category check
            if (
                not (
                    self.pedigree[sample_id].affected == PEDDY_AFFECTED
                    and principal.sample_category_check(sample_id, allow_support=False)
                )
            ) or (
                principal.depths[sample_id] < self.minimum_depth
                and not principal.info.get('categoryboolean1')
            ):
                continue

            # check if this is a possible candidate for homozygous inheritance
            if not self.check_familial_inheritance(
                sample_id=sample_id,
                called_variants=males,
                partial_pen=partial_pen,
            ):
                continue

            classifications.append(
                ReportedVariant(
                    sample=sample_id,
                    family=self.pedigree[sample_id].family_id,
                    gene=principal.info.get('gene_id'),
                    var_data=MinimalVariant(principal, sample_id),
                    genotypes=self.get_family_genotypes(
                        variant=principal, sample_id=sample_id
                    ),
                    reasons={self.applied_moi},
                    flags=principal.get_sample_flags(sample_id),
                )
            )
        return classifications


class XRecessiveFemaleHom(BaseMoi):
    """
    only consider HOM females
    """

    def __init__(
        self,
        pedigree: Ped,
        applied_moi: str = 'X_Recessive HOM Female',
    ):
        """
        set parameters specific to recessive tests

        Args:
            pedigree ():
            applied_moi ():
        """

        self.hom_rec_threshold = get_config()['moi_tests'][GNOMAD_REC_HOM_THRESHOLD]
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: AbstractVariant,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,
    ) -> list[ReportedVariant]:
        """

        Args:
            principal ():
            comp_het ():
            partial_pen ():
        """

        classifications = []

        # remove from analysis if too many homs are present in population databases
        if principal.support_only or (
            (
                self.check_frequency(
                    principal.info, threshold=self.hom_rec_threshold, keys=INFO_HOMS
                )
            )
            and not principal.info.get('categoryboolean1')
        ):
            return classifications

        # never consider support homs
        samples_to_check = {
            sam for sam in principal.hom_samples if self.pedigree[sam].sex == 'female'
        }

        for sample_id in samples_to_check:

            # specific affected sample category check
            if (
                not (
                    self.pedigree[sample_id].affected == PEDDY_AFFECTED
                    and principal.sample_category_check(sample_id, allow_support=False)
                )
            ) or (
                principal.depths[sample_id] < self.minimum_depth
                and not principal.info.get('categoryboolean1')
            ):
                continue

            # check if this is a possible candidate for homozygous inheritance
            if not self.check_familial_inheritance(
                sample_id=sample_id,
                called_variants=samples_to_check,
                partial_pen=partial_pen,
            ):
                continue

            classifications.append(
                ReportedVariant(
                    sample=sample_id,
                    family=self.pedigree[sample_id].family_id,
                    gene=principal.info.get('gene_id'),
                    var_data=MinimalVariant(principal, sample_id),
                    genotypes=self.get_family_genotypes(
                        variant=principal, sample_id=sample_id
                    ),
                    reasons={self.applied_moi},
                    flags=principal.get_sample_flags(sample_id),
                )
            )
        return classifications


class XRecessiveFemaleCH(BaseMoi):
    """
    ignore males, accept female comp-het only
    """

    def __init__(
        self,
        pedigree: Ped,
        applied_moi: str = 'X_RecessiveFemaleCompHet',
    ):
        """
        set parameters specific to recessive tests

        Args:
            pedigree ():
            applied_moi ():
        """

        self.hom_rec_threshold = get_config()['moi_tests'][GNOMAD_REC_HOM_THRESHOLD]
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: AbstractVariant,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,
    ) -> list[ReportedVariant]:
        """

        Args:
            principal ():
            comp_het ():
            partial_pen ():
        """

        if comp_het is None:
            comp_het = {}

        classifications = []

        # remove from analysis if too many homs are present in population databases
        if self.check_frequency(
            principal.info, threshold=self.hom_rec_threshold, keys=INFO_HOMS
        ) and not principal.info.get('categoryboolean1'):
            return classifications

        het_females = {
            sam for sam in principal.het_samples if self.pedigree[sam].sex == 'female'
        }

        # if het females are present, try and find support
        for sample_id in het_females:

            # don't run primary analysis for unaffected
            # we require this specific sample to be categorised - check Cat 4 contents
            if (
                not (
                    self.pedigree[sample_id].affected == PEDDY_AFFECTED
                    and principal.sample_category_check(sample_id, allow_support=True)
                )
            ) or (
                principal.depths[sample_id] < self.minimum_depth
                and not principal.info.get('categoryboolean1')
            ):
                continue

            for partner in check_for_second_hit(
                first_variant=principal.coords.string_format,
                comp_hets=comp_het,
                sample=sample_id,
            ):

                # allow for de novo check - also screen out high-AF partners
                if (
                    not partner.sample_category_check(sample_id, allow_support=True)
                    or (
                        self.check_frequency(
                            partner.info,
                            threshold=self.hom_rec_threshold,
                            keys=INFO_HOMS,
                        )
                        and not partner.info.get('categoryboolean1')
                    )
                ) or (
                    principal.sample_support_only(sample_id)
                    and partner.sample_support_only(sample_id)
                ):
                    continue

                # check for minimum depth in partner
                if partner.depths[
                    sample_id
                ] < self.minimum_depth and not partner.info.get('categoryboolean1'):
                    continue

                if not self.check_comp_het(
                    sample_id=sample_id, variant_1=principal, variant_2=partner
                ):
                    continue

                classifications.append(
                    ReportedVariant(
                        sample=sample_id,
                        family=self.pedigree[sample_id].family_id,
                        gene=principal.info.get('gene_id'),
                        var_data=MinimalVariant(principal, sample_id),
                        reasons={self.applied_moi},
                        genotypes=self.get_family_genotypes(
                            variant=principal, sample_id=sample_id
                        ),
                        support_vars={partner.coords.string_format},
                        flags=principal.get_sample_flags(sample_id)
                        + partner.get_sample_flags(sample_id),
                    )
                )

        return classifications
