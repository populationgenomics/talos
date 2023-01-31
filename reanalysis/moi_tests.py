"""
A number of classes, each representing one Mode of Inheritance
One class (MoiRunner) to run all the appropriate MOIs on a variant

Reduce the PanelApp plain text MOI description into a few categories
We then run a permissive MOI match for the variant,

This model does not apply anything to MT, I expect those to default
to a Monogenic MOI
"""


from abc import abstractmethod
from copy import deepcopy

from peddy.peddy import Ped, PHENOTYPE

from cpg_utils.config import get_config

from reanalysis.utils import AbstractVariant, CompHetDict, ReportedVariant, X_CHROMOSOME


# config keys to use for dominant MOI tests
GNOMAD_RARE_THRESHOLD = 'gnomad_dominant'
GNOMAD_AD_AC_THRESHOLD = 'gnomad_max_ac_dominant'
GNOMAD_DOM_HOM_THRESHOLD = 'gnomad_max_homs_dominant'
GNOMAD_REC_HOM_THRESHOLD = 'gnomad_max_homs_recessive'
GNOMAD_HEMI_THRESHOLD = 'gnomad_max_hemi'
INFO_HOMS = {'gnomad_hom', 'gnomad_ex_hom'}
INFO_HEMI = {'gnomad_hemi', 'gnomad_ex_hemi'}
PEDDY_AFFECTED = PHENOTYPE().AFFECTED
VAR_FIELDS_TO_REMOVE = [
    'het_samples',
    'hom_samples',
    'boolean_categories',
    'sample_categories',
    'sample_support',
    'ab_ratios',
]


def minimise_variant(variant: AbstractVariant, sample_id: str) -> AbstractVariant:
    """
    When we find a matching MOI, store a sample-specific
    duplicate of the variant details
    Args:
        variant ():
        sample_id ():

    Returns:

    """
    var_copy = deepcopy(variant)
    var_copy.categories = var_copy.category_values(sample=sample_id)
    for key in VAR_FIELDS_TO_REMOVE:
        if key in var_copy.__dict__:
            del var_copy.__dict__[key]
    return var_copy


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
                RecessiveAutosomal(pedigree=pedigree),
            ]
        elif target_moi == 'Biallelic':
            self.filter_list = [RecessiveAutosomal(pedigree=pedigree)]

        elif target_moi == 'Hemi_Mono_In_Female':
            self.filter_list = [
                XRecessive(pedigree=pedigree),
                XDominant(pedigree=pedigree),
            ]

        elif target_moi == 'Hemi_Bi_In_Female':
            self.filter_list = [XRecessive(pedigree=pedigree)]

        else:
            raise Exception(f'MOI type {target_moi} is not addressed in MOI')

    def run(
        self,
        principal_var,
        comp_het: CompHetDict | None = None,
        partial_penetrance: bool = False,
    ) -> list[ReportedVariant]:
        """
        run method - triggers each relevant inheritance model
        :param principal_var: the variant we are focused on
        :param comp_het:
        :param partial_penetrance:
        :return:
        """

        if comp_het is None:
            comp_het = {}

        moi_matched = []
        for model in self.filter_list:
            moi_matched.extend(
                model.run(
                    principal_var=principal_var,
                    comp_het=comp_het,
                    partial_penetrance=partial_penetrance,
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
            raise Exception('An applied MOI needs to reach the Base Class')
        self.pedigree = pedigree
        self.applied_moi = applied_moi

    @abstractmethod
    def run(
        self,
        principal_var: AbstractVariant,
        comp_het: CompHetDict | None = None,
        partial_penetrance: bool = False,
    ) -> list[ReportedVariant]:
        """
        run all applicable inheritance patterns and finds good fits
        :param principal_var:
        :param comp_het: dictionary of compound-hets
        :param partial_penetrance: default to False, permit setting True
        :return:
        """

    def check_familial_inheritance(
        self, sample_id: str, called_variants: set[str], partial_pen: bool = False
    ) -> bool:
        """
        sex-agnostic check for single variant inheritance
        designed to be called recursively, this method will take a single sample
        as a starting point and check that the MOI is viable across the whole family
        *we are assuming complete penetrance with this version*

        - find the family ID from this sample
        - iterate through all family members, and check that MOI holds for all

        Return False if a participant fails tests, True if all checked are ok

        :param sample_id:
        :param called_variants: the set of sample_ids which have this variant
        :param partial_pen: if True, permit unaffected has variant call

        NOTE: this called_variants pool is prepared before calling this method.
        If we only want to check hom calls, only send hom calls. If we are checking for
        dominant conditions, we would bundle hom and het calls into this set

        This check should work for both autosomal & sex chromosomes...
        for X-chrom biallelic, send female Homs, and all male variants
        For autosomal biallelic send male and female Homs
        :return:
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

    def check_familial_comp_het(
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

        :param pedigree: not yet implemented
        :param applied_moi:
        """
        self.ad_threshold = get_config()['moi_tests'][GNOMAD_RARE_THRESHOLD]
        self.ac_threshold = get_config()['moi_tests'][GNOMAD_AD_AC_THRESHOLD]
        self.hom_threshold = get_config()['moi_tests'][GNOMAD_DOM_HOM_THRESHOLD]
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal_var: AbstractVariant,
        comp_het: CompHetDict | None = None,
        partial_penetrance: bool = False,
    ) -> list[ReportedVariant]:
        """
        Simplest MOI, exclusions based on HOM count and AF
        Args:
            principal_var ():
            comp_het ():
            partial_penetrance ():
        """

        classifications = []

        # more stringent Pop.Freq checks for dominant
        if (
            principal_var.info.get('gnomad_af', 0) > self.ad_threshold
            or any(
                {
                    principal_var.info.get(hom_key, 0) > self.hom_threshold
                    for hom_key in INFO_HOMS
                }
            )
            or principal_var.info.get('gnomad_ac', 0) > self.ac_threshold
        ) and not principal_var.info.get('categoryboolean1'):
            return classifications

        # autosomal dominant doesn't require support, but consider het and hom
        samples_with_this_variant = principal_var.het_samples.union(
            principal_var.hom_samples
        )
        for sample_id in samples_with_this_variant:

            # skip primary analysis for unaffected members
            if not self.pedigree[sample_id].affected == PEDDY_AFFECTED:
                continue

            # we require this specific sample to be categorised - check Cat 4 contents
            if not principal_var.sample_specific_category_check(sample_id):
                continue

            # check if this is a candidate for dominant inheritance
            if not self.check_familial_inheritance(
                sample_id=sample_id,
                called_variants=samples_with_this_variant,
                partial_pen=partial_penetrance,
            ):
                continue

            var_copy = minimise_variant(variant=principal_var, sample_id=sample_id)
            classifications.append(
                ReportedVariant(
                    sample=sample_id,
                    family=self.pedigree[sample_id].family_id,
                    gene=var_copy.info.get('gene_id'),
                    var_data=var_copy,
                    reasons={self.applied_moi},
                    genotypes=self.get_family_genotypes(
                        variant=principal_var, sample_id=sample_id
                    ),
                    flags=principal_var.get_sample_flags(sample_id),
                )
            )

        return classifications


class RecessiveAutosomal(BaseMoi):
    """
    inheritance test class for Recessive inheritance
    requires single hom variant, or compound het
    """

    def __init__(
        self,
        pedigree: Ped,
        applied_moi: str = 'Autosomal Recessive',
    ):
        """ """
        self.hom_threshold = get_config()['moi_tests'][GNOMAD_REC_HOM_THRESHOLD]
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal_var: AbstractVariant,
        comp_het: CompHetDict | None = None,
        partial_penetrance: bool = False,
    ) -> list[ReportedVariant]:
        """
        valid if present as hom, or compound het
        counts as being phased if a compound het is split between parents
        Clarify if we want to consider a homozygous variant as 2 hets

        Args:
            principal_var (AbstractVariant): main variant being evaluated
            comp_het (dict): comp-het partners
            partial_penetrance (bool):

        Returns:
            list[ReportedVariant]: data object if RecessiveAutosomal fits
        """

        if comp_het is None:
            comp_het = {}

        classifications = []

        # remove if too many homs are present in population databases
        # no stricter AF here - if we choose to, we can apply while labelling
        if any(
            {
                principal_var.info.get(hom_key, 0) > self.hom_threshold
                for hom_key in INFO_HOMS
            }
        ) and not principal_var.info.get('categoryboolean1'):
            return classifications

        # homozygous is relevant directly
        for sample_id in principal_var.hom_samples:

            # skip primary analysis for unaffected members
            if not self.pedigree[sample_id].affected == PEDDY_AFFECTED:
                continue

            # we require this specific sample to be categorised - check Cat 4 contents
            # this shouldn't be possible, de novo homozygous?!
            if not principal_var.sample_specific_category_check(sample_id):
                continue

            # check if this is a possible candidate for homozygous inheritance
            if not self.check_familial_inheritance(
                sample_id=sample_id,
                called_variants=principal_var.hom_samples,
                partial_pen=partial_penetrance,
            ):
                continue

            var_copy = minimise_variant(variant=principal_var, sample_id=sample_id)
            classifications.append(
                ReportedVariant(
                    sample=sample_id,
                    family=self.pedigree[sample_id].family_id,
                    gene=var_copy.info.get('gene_id'),
                    var_data=var_copy,
                    genotypes=self.get_family_genotypes(
                        variant=principal_var, sample_id=sample_id
                    ),
                    reasons={f'{self.applied_moi} Homozygous'},
                    flags=principal_var.get_sample_flags(sample_id),
                )
            )

        # if hets are present, try and find support
        for sample_id in principal_var.het_samples:

            # skip primary analysis for unaffected members
            if not self.pedigree[sample_id].affected == PEDDY_AFFECTED:
                continue

            # we require this specific sample to be categorised - check Cat 4 contents
            if not principal_var.sample_specific_category_check(sample_id):
                continue

            for partner_variant in check_for_second_hit(
                first_variant=principal_var.coords.string_format,
                comp_hets=comp_het,
                sample=sample_id,
            ):

                # categorised for this specific sample, allow support in partner
                # - also screen out high-AF partners
                if (
                    not (
                        partner_variant.sample_specific_category_check(sample_id)
                        or partner_variant.has_support
                    )
                    or any(  # allow for clinvar here
                        {
                            partner_variant.info.get(hom_key, 0) > self.hom_threshold
                            for hom_key in INFO_HOMS
                        }
                    )
                    and not principal_var.info.get('categoryboolean1')
                ):
                    continue

                # check if this is a candidate for comp-het inheritance
                if not self.check_familial_comp_het(
                    sample_id=sample_id,
                    variant_1=principal_var,
                    variant_2=partner_variant,
                ):
                    continue

                var_copy = minimise_variant(variant=principal_var, sample_id=sample_id)
                classifications.append(
                    ReportedVariant(
                        sample=sample_id,
                        family=self.pedigree[sample_id].family_id,
                        gene=var_copy.info.get('gene_id'),
                        var_data=var_copy,
                        reasons={f'{self.applied_moi} Compound-Het'},
                        supported=True,
                        genotypes=self.get_family_genotypes(
                            variant=principal_var, sample_id=sample_id
                        ),
                        support_vars=[partner_variant.coords.string_format],
                        flags=principal_var.get_sample_flags(sample_id)
                        + partner_variant.get_sample_flags(sample_id),
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
        principal_var: AbstractVariant,
        comp_het: CompHetDict | None = None,
        partial_penetrance: bool = False,
    ) -> list[ReportedVariant]:
        """
        if variant is present and sufficiently rare, we take it

        :param principal_var:
        :param comp_het:
        :param partial_penetrance:
        :return:
        """

        classifications = []

        if principal_var.coords.chrom.lower() != 'x':
            raise Exception(
                f'X-Chromosome MOI given for variant on {principal_var.coords.chrom}'
            )

        # more stringent Pop.Freq checks for dominant - hemi restriction
        if (
            principal_var.info.get('gnomad_af', 0) > self.ad_threshold
            or any(
                {
                    principal_var.info.get(hom_key, 0) > self.hom_threshold
                    for hom_key in INFO_HOMS
                }
            )
            or principal_var.info.get('gnomad_ac', 0) > self.ac_threshold
            or any(
                {
                    principal_var.info.get(hemi_key, 0) > self.hemi_threshold
                    for hemi_key in INFO_HEMI
                }
            )
        ) and not principal_var.info.get('categoryboolean1'):
            return classifications

        # all samples which have a variant call
        samples_with_this_variant = principal_var.het_samples.union(
            principal_var.hom_samples
        )

        for sample_id in samples_with_this_variant:

            # skip primary analysis for unaffected members
            if not self.pedigree[sample_id].affected == PEDDY_AFFECTED:
                continue

            # we require this specific sample to be categorised - check Cat 4 contents
            if not principal_var.sample_specific_category_check(sample_id):
                continue

            # check if this is a candidate for dominant inheritance
            if not self.check_familial_inheritance(
                sample_id=sample_id,
                called_variants=samples_with_this_variant,
                partial_pen=partial_penetrance,
            ):
                continue

            # passed inheritance test, create the record

            var_copy = minimise_variant(variant=principal_var, sample_id=sample_id)
            classifications.append(
                ReportedVariant(
                    sample=sample_id,
                    family=self.pedigree[sample_id].family_id,
                    gene=var_copy.info.get('gene_id'),
                    var_data=var_copy,
                    reasons={
                        f'{self.applied_moi} '
                        f'{self.pedigree[sample_id].sex.capitalize()}'
                    },
                    genotypes=self.get_family_genotypes(
                        variant=principal_var, sample_id=sample_id
                    ),
                    flags=principal_var.get_sample_flags(sample_id),
                )
            )
        return classifications


class XRecessive(BaseMoi):
    """
    for males accept het** - male variants HOM because GATK
    effectively the same as AutosomalDominant?
    """

    def __init__(
        self,
        pedigree: Ped,
        applied_moi: str = 'X_Recessive',
    ):
        """
        set parameters specific to recessive tests
        and take the comp-het dictionary
        :param pedigree:
        :param applied_moi:
        """

        self.hom_dom_threshold = get_config()['moi_tests'][GNOMAD_DOM_HOM_THRESHOLD]
        self.hom_rec_threshold = get_config()['moi_tests'][GNOMAD_REC_HOM_THRESHOLD]
        self.hemi_threshold = get_config()['moi_tests'][GNOMAD_HEMI_THRESHOLD]

        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal_var: AbstractVariant,
        comp_het: CompHetDict | None = None,
        partial_penetrance: bool = False,
    ) -> list[ReportedVariant]:
        """

        :param principal_var:
        :param comp_het:
        :param partial_penetrance:
        :return:
        """

        if comp_het is None:
            comp_het = {}

        classifications = []

        # remove from analysis if too many homs are present in population databases
        if any(
            {
                principal_var.info.get(hom_key, 0) > self.hom_dom_threshold
                for hom_key in INFO_HOMS
            }
        ) and not principal_var.info.get('categoryboolean1'):
            return classifications

        # X-relevant, we separate male & females
        # combine het and hom here, we don't trust the variant callers
        # if hemi count is too high, don't consider males
        males = {
            sam
            for sam in principal_var.het_samples.union(principal_var.hom_samples)
            if self.pedigree[sam].sex == 'male'
        }

        # split female calls into 2 categories
        het_females = {
            sam
            for sam in principal_var.het_samples
            if self.pedigree[sam].sex == 'female'
        }
        hom_females = {
            sam
            for sam in principal_var.hom_samples
            if self.pedigree[sam].sex == 'female'
        }

        # if het females are present, try and find support
        for sample_id in het_females:

            # don't run primary analysis for unaffected
            # we require this specific sample to be categorised - check Cat 4 contents
            if not (
                self.pedigree[sample_id].affected == PEDDY_AFFECTED
                and principal_var.sample_specific_category_check(sample_id)
            ):
                continue

            for partner_variant in check_for_second_hit(
                first_variant=principal_var.coords.string_format,
                comp_hets=comp_het,
                sample=sample_id,
            ):

                # allow for de novo check - also screen out high-AF partners
                if (
                    not (
                        partner_variant.sample_specific_category_check(sample_id)
                        or partner_variant.has_support
                    )
                    or any(
                        {
                            partner_variant.info.get(hom_key, 0)
                            >= self.hom_rec_threshold
                            for hom_key in INFO_HOMS
                        }
                    )
                    and not principal_var.info.get('categoryboolean1')
                ):
                    continue

                if not self.check_familial_comp_het(
                    sample_id=sample_id,
                    variant_1=principal_var,
                    variant_2=partner_variant,
                ):
                    continue

                var_copy = minimise_variant(variant=principal_var, sample_id=sample_id)
                classifications.append(
                    ReportedVariant(
                        sample=sample_id,
                        family=self.pedigree[sample_id].family_id,
                        gene=var_copy.info.get('gene_id'),
                        var_data=var_copy,
                        reasons={f'{self.applied_moi} Compound-Het Female'},
                        supported=True,
                        genotypes=self.get_family_genotypes(
                            variant=principal_var, sample_id=sample_id
                        ),
                        support_vars=[partner_variant.coords.string_format],
                        flags=principal_var.get_sample_flags(sample_id)
                        + partner_variant.get_sample_flags(sample_id),
                    )
                )

        # remove from analysis if too many homs are present in population databases
        if any(
            {
                principal_var.info.get(hom_key, 0) > self.hom_rec_threshold
                for hom_key in INFO_HOMS
            }
        ) and not principal_var.info.get('categoryboolean1'):
            return classifications

        # find all het males and hom females
        samples_to_check = males.union(hom_females)
        for sample_id in samples_to_check:

            # specific affected sample category check
            if not (
                self.pedigree[sample_id].affected == PEDDY_AFFECTED
                and principal_var.sample_specific_category_check(sample_id)
            ):
                continue

            # if this is male, and hemi count is high, skip the sample
            if self.pedigree[sample_id].sex == 'male' and any(
                {
                    principal_var.info.get(hemi_key, 0) > self.hemi_threshold
                    for hemi_key in INFO_HEMI
                }
            ):
                continue

            # check if this is a possible candidate for homozygous inheritance
            if not self.check_familial_inheritance(
                sample_id=sample_id,
                called_variants=samples_to_check,
                partial_pen=partial_penetrance,
            ):
                continue

            var_copy = minimise_variant(variant=principal_var, sample_id=sample_id)
            classifications.append(
                ReportedVariant(
                    sample=sample_id,
                    family=self.pedigree[sample_id].family_id,
                    gene=var_copy.info.get('gene_id'),
                    var_data=var_copy,
                    genotypes=self.get_family_genotypes(
                        variant=principal_var, sample_id=sample_id
                    ),
                    reasons={
                        f'{self.applied_moi} '
                        f'{self.pedigree[sample_id].sex.capitalize()}'
                    },
                    flags=principal_var.get_sample_flags(sample_id),
                )
            )
        return classifications
