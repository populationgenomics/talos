"""
A number of classes, each representing one Mode of Inheritance
One class (MoiRunner) to run all the appropriate MOIs on a variant

Reduce the PanelApp plain text MOI description into a few categories
We then run a permissive MOI match for the variant
"""
# mypy: ignore-errors
from abc import abstractmethod

from cpg_utils.config import get_config
from peddy.peddy import Ped, PHENOTYPE

from reanalysis.models import (
    SmallVariant,
    StructuralVariant,
    ReportVariant,
    VARIANT_MODELS,
)
from reanalysis.utils import CompHetDict, X_CHROMOSOME

# config keys to use for dominant MOI tests
CALLSET_AF_SV_DOMINANT = 'callset_af_sv_dominant'
GNOMAD_RARE_THRESHOLD = 'gnomad_dominant'
GNOMAD_AD_AC_THRESHOLD = 'gnomad_max_ac_dominant'
GNOMAD_DOM_HOM_THRESHOLD = 'gnomad_max_homs_dominant'
GNOMAD_REC_HOM_THRESHOLD = 'gnomad_max_homs_recessive'
GNOMAD_HEMI_THRESHOLD = 'gnomad_max_hemi'
INFO_HOMS = {'gnomad_hom', 'gnomad_ex_hom'}
INFO_HEMI = {'gnomad_hemi', 'gnomad_ex_hemi'}
PEDDY_AFFECTED = PHENOTYPE().AFFECTED
SV_AF_KEY = 'gnomad_v2.1_sv_AF'
SV_HEMI = {'male_n_hemialt'}
SV_HOMS = {'male_n_homalt', 'female_n_homalt'}


def check_for_second_hit(
    first_variant: str,
    comp_hets: CompHetDict,
    sample: str,
    require_non_support: bool = False,
) -> list[VARIANT_MODELS]:
    """
    checks for a second hit partner in this gene

    Example formatting of the comp-het dict
    {
        "SampleID": {
            "12-52287177-T-C": [
                Variant(12-52287180-TGG-T)
            ],
            "12-52287180-TGG-T": [
                Variant(12-52287177-T-C)
            ]
        } ...
    }

    Args:
        first_variant (str): string representation of variant1
        comp_hets (dict[str, Variant]): lookup for compound hets
        sample (str): sample ID
        require_non_support (bool): if true, don't return Support only

    Returns:
        a list of variants which are potential partners
    """

    # check if the sample has any comp-hets
    if sample not in comp_hets:
        return []

    partners = comp_hets[sample].get(first_variant, [])
    if require_non_support:
        return [
            partner for partner in partners if not partner.sample_support_only(sample)
        ]
    else:
        return partners


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
    ) -> list[ReportVariant]:
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
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,
    ) -> list[ReportVariant]:
        """
        run all applicable inheritance patterns and finds good fits
        """

    def check_affected_category_depths(
        self, variant: VARIANT_MODELS, sample_id: str
    ) -> bool:
        """

        Args:
            variant ():
            sample_id ():

        Returns:

        """
        if (
            not (
                self.pedigree[sample_id].affected == PEDDY_AFFECTED
                and variant.sample_category_check(sample_id, allow_support=False)
            )
        ) or variant.check_read_depth(
            sample_id,
            self.minimum_depth,
            var_is_cat_1=variant.info.get('categoryboolean1'),
        ):
            return True
        return False

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
        self, variant: VARIANT_MODELS, sample_id: str
    ) -> dict[str, str]:
        """

        Args:
            variant (VARIANT_MODELS):
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

            if variant.coordinates.chrom in X_CHROMOSOME:
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
    def check_frequency_passes(info: dict, thresholds: dict[str, int | float]) -> bool:
        """
        Method to check multiple info keys against a single threshold
        This just reduces the line count, as this is called a bunch of times

        Args:
            info (): the dict of values for this dict
            thresholds (): the dict of keys - thresholds to test against

        Returns:
            True if any of the info attributes is above the threshold
        """
        return all({info.get(key, 0) <= test for key, test in thresholds.items()})

    def check_comp_het(
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
        self.hom_threshold = get_config()['moi_tests'][GNOMAD_DOM_HOM_THRESHOLD]
        self.sv_af_threshold = get_config()['moi_tests'][CALLSET_AF_SV_DOMINANT]

        # prepare the AF test dicts
        self.freq_tests = {
            SmallVariant.__name__: {key: self.hom_threshold for key in INFO_HOMS}
            | {
                'gnomad_ac': self.ac_threshold,
                'gnomad_af': self.ad_threshold,
            },
            StructuralVariant.__name__: {
                'af': self.sv_af_threshold,
                SV_AF_KEY: self.sv_af_threshold,
            },
        }
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,
    ) -> list[ReportVariant]:
        """
        Simplest MOI, exclusions based on HOM count and AF
        Args:
            principal ():
            comp_het ():
            partial_pen ():
        """

        classifications = []

        # reject support for dominant MOI, apply checks based on var type
        if principal.support_only or not (
            self.check_frequency_passes(
                principal.info, self.freq_tests[principal.__class__.__name__]
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
                principal.check_read_depth(
                    sample_id,
                    self.minimum_depth,
                    var_is_cat_1=principal.info.get('categoryboolean1'),
                )
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
                ReportVariant(
                    sample=sample_id,
                    family=self.pedigree[sample_id].family_id,
                    gene=principal.info.get('gene_id'),
                    var_data=principal,
                    categories=principal.category_values(sample_id),
                    reasons={self.applied_moi},
                    genotypes=self.get_family_genotypes(
                        variant=principal, sample_id=sample_id
                    ),
                    flags=principal.get_sample_flags(sample_id),
                    independent=True,
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
        self.hom_threshold = get_config()['moi_tests'][GNOMAD_REC_HOM_THRESHOLD]
        self.freq_tests = {
            SmallVariant.__name__: {key: self.hom_threshold for key in INFO_HOMS},
            StructuralVariant.__name__: {key: self.hom_threshold for key in SV_HOMS},
        }
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,
    ) -> list[ReportVariant]:
        """
        valid if present as compound het
        counts as being phased if a compound het is split between parents

        Args:
            principal (VARIANT_MODELS): main variant being evaluated
            comp_het (dict): comp-het partners
            partial_pen (bool):

        Returns:
            list[ReportVariant]: data object if RecessiveAutosomal fits
        """

        if comp_het is None:
            comp_het = {}

        classifications = []

        # remove if too many homs are present in population databases
        if not (
            self.check_frequency_passes(
                principal.info, self.freq_tests[principal.__class__.__name__]
            )
            or principal.info.get('categoryboolean1')
        ):
            return classifications

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
                principal.check_read_depth(
                    sample_id,
                    self.minimum_depth,
                    principal.info.get('categoryboolean1'),
                )
            ):
                continue

            for partner_variant in check_for_second_hit(
                first_variant=principal.coordinates.string_format,
                comp_hets=comp_het,
                sample=sample_id,
                require_non_support=principal.sample_support_only(sample_id),
            ):

                if partner_variant.check_read_depth(
                    sample_id,
                    self.minimum_depth,
                    partner_variant.info.get('categoryboolean1'),
                ) or not (
                    self.check_frequency_passes(
                        partner_variant.info,
                        self.freq_tests[partner_variant.__class__.__name__],
                    )
                ):
                    continue

                # # categorised for this specific sample, allow support in partner
                # # - also screen out high-AF partners
                # if not partner_variant.sample_category_check(
                #     sample_id, allow_support=True
                # ):
                #     continue

                # check if this is a candidate for comp-het inheritance
                if not self.check_comp_het(
                    sample_id=sample_id, variant_1=principal, variant_2=partner_variant
                ):
                    continue

                classifications.append(
                    ReportVariant(
                        sample=sample_id,
                        family=self.pedigree[sample_id].family_id,
                        gene=principal.info.get('gene_id'),
                        var_data=principal,
                        categories=principal.category_values(sample_id),
                        reasons={self.applied_moi},
                        genotypes=self.get_family_genotypes(
                            variant=principal, sample_id=sample_id
                        ),
                        support_vars={partner_variant.coordinates.string_format},
                        flags=principal.get_sample_flags(sample_id)
                        | partner_variant.get_sample_flags(sample_id),
                        independent=False,
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
        self.freq_tests = {
            SmallVariant.__name__: {key: self.hom_threshold for key in INFO_HOMS},
            StructuralVariant.__name__: {key: self.hom_threshold for key in SV_HOMS},
        }
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,
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
        if principal.support_only or not (
            self.check_frequency_passes(
                principal.info, self.freq_tests[principal.__class__.__name__]
            )
            or principal.info.get('categoryboolean1')
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
            ) or principal.check_read_depth(
                sample_id, self.minimum_depth, principal.info.get('categoryboolean1')
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
                ReportVariant(
                    sample=sample_id,
                    family=self.pedigree[sample_id].family_id,
                    gene=principal.info.get('gene_id'),
                    var_data=principal,
                    categories=principal.category_values(sample_id),
                    genotypes=self.get_family_genotypes(
                        variant=principal, sample_id=sample_id
                    ),
                    reasons={self.applied_moi},
                    flags=principal.get_sample_flags(sample_id),
                    independent=True,
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

        self.freq_tests = {
            SmallVariant.__name__: {key: self.hom_threshold for key in INFO_HOMS}
            | {key: self.hemi_threshold for key in INFO_HEMI}
            | {
                'gnomad_ac': self.ac_threshold,
                'gnomad_af': self.ad_threshold,
            },
            StructuralVariant.__name__: {key: self.hom_threshold for key in SV_HOMS}
            | {key: self.hemi_threshold for key in SV_HEMI},
        }

        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,
    ) -> list[ReportVariant]:
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
        if not self.check_frequency_passes(
            principal.info, self.freq_tests[principal.__class__.__name__]
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
            ) or principal.check_read_depth(
                sample_id, self.minimum_depth, principal.info.get('categoryboolean1')
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
                ReportVariant(
                    sample=sample_id,
                    family=self.pedigree[sample_id].family_id,
                    gene=principal.info.get('gene_id'),
                    var_data=principal,
                    categories=principal.category_values(sample_id),
                    reasons={self.applied_moi},
                    genotypes=self.get_family_genotypes(
                        variant=principal, sample_id=sample_id
                    ),
                    flags=principal.get_sample_flags(sample_id),
                    independent=True,
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

        self.freq_tests = {
            SmallVariant.__name__: {key: self.hom_dom_threshold for key in INFO_HOMS}
            | {key: self.hemi_threshold for key in INFO_HEMI},
            StructuralVariant.__name__: {key: self.hom_dom_threshold for key in SV_HOMS}
            | {key: self.hemi_threshold for key in SV_HEMI},
        }

        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,
    ) -> list[ReportVariant]:
        """
        Args:
            principal ():
            comp_het ():
            partial_pen ():
        """

        classifications = []

        # remove from analysis if too many homs are present in population databases
        if not self.check_frequency_passes(
            principal.info, self.freq_tests[principal.__class__.__name__]
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

            # specific affected sample category check, never consider support on X for males
            if (
                not (
                    self.pedigree[sample_id].affected == PEDDY_AFFECTED
                    and principal.sample_category_check(sample_id, allow_support=False)
                )
            ) or principal.check_read_depth(
                sample_id, self.minimum_depth, principal.info.get('categoryboolean1')
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
                ReportVariant(
                    sample=sample_id,
                    family=self.pedigree[sample_id].family_id,
                    gene=principal.info.get('gene_id'),
                    var_data=principal,
                    categories=principal.category_values(sample_id),
                    genotypes=self.get_family_genotypes(
                        variant=principal, sample_id=sample_id
                    ),
                    reasons={self.applied_moi},
                    flags=principal.get_sample_flags(sample_id),
                    independent=True,
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
        self.freq_tests = {
            SmallVariant.__name__: {key: self.hom_rec_threshold for key in INFO_HOMS},
            StructuralVariant.__name__: {
                key: self.hom_rec_threshold for key in SV_HOMS
            },
        }
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,
    ) -> list[ReportVariant]:
        """

        Args:
            principal ():
            comp_het ():
            partial_pen ():
        """

        classifications = []

        # remove from analysis if too many homs are present in population databases
        if principal.support_only or not (
            self.check_frequency_passes(
                principal.info, self.freq_tests[principal.__class__.__name__]
            )
            or principal.info.get('categoryboolean1')
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
            ) or principal.check_read_depth(
                sample_id, self.minimum_depth, principal.info.get('categoryboolean1')
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
                ReportVariant(
                    sample=sample_id,
                    family=self.pedigree[sample_id].family_id,
                    gene=principal.info.get('gene_id'),
                    var_data=principal,
                    categories=principal.category_values(sample_id),
                    genotypes=self.get_family_genotypes(
                        variant=principal, sample_id=sample_id
                    ),
                    reasons={self.applied_moi},
                    flags=principal.get_sample_flags(sample_id),
                    independent=True,
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
        self.freq_tests = {
            SmallVariant.__name__: {key: self.hom_rec_threshold for key in INFO_HOMS},
            StructuralVariant.__name__: {
                key: self.hom_rec_threshold for key in SV_HOMS
            },
        }
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,
    ) -> list[ReportVariant]:
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
        if not (
            self.check_frequency_passes(
                principal.info, self.freq_tests[principal.__class__.__name__]
            )
            or principal.info.get('categoryboolean1')
        ):
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
            ) or principal.check_read_depth(
                sample_id, self.minimum_depth, principal.info.get('categoryboolean1')
            ):
                continue

            for partner in check_for_second_hit(
                first_variant=principal.coordinates.string_format,
                comp_hets=comp_het,
                sample=sample_id,
                require_non_support=principal.sample_support_only(sample_id),
            ):

                # allow for de novo check - also screen out high-AF partners
                if not partner.sample_category_check(
                    sample_id, allow_support=True
                ) or not (
                    self.check_frequency_passes(
                        partner.info, self.freq_tests[partner.__class__.__name__]
                    )
                    or partner.info.get('categoryboolean1')
                ):
                    continue

                # check for minimum depth in partner
                if partner.check_read_depth(
                    sample_id, self.minimum_depth, partner.info.get('categoryboolean1')
                ):
                    continue

                if not self.check_comp_het(
                    sample_id=sample_id, variant_1=principal, variant_2=partner
                ):
                    continue

                classifications.append(
                    ReportVariant(
                        sample=sample_id,
                        family=self.pedigree[sample_id].family_id,
                        gene=principal.info.get('gene_id'),
                        var_data=principal,
                        categories=principal.category_values(sample_id),
                        reasons={self.applied_moi},
                        genotypes=self.get_family_genotypes(
                            variant=principal, sample_id=sample_id
                        ),
                        support_vars={partner.coordinates.string_format},
                        flags=principal.get_sample_flags(sample_id)
                        | partner.get_sample_flags(sample_id),
                        independent=False,
                    )
                )

        return classifications
