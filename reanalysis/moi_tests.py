"""
A number of classes, each representing one Mode of Inheritance
One class (MoiRunner) to run all the appropriate MOIs on a variant

Reduce the PanelApp plain text MOI description into a few categories
We then run a permissive MOI match for the variant,
e.g. if the MOI is Dominant, we may also be interested in Recessive (?)
e.g. if the MOI is X-linked Dom, we also search X-linked Recessive (?)
    - relevant for females, if the pedigree is loaded
    - also depends on the accuracy of a male hemi call

This will come down to clinician preference, some may exclusively
prefer dominant model = monogenic search

This model does not apply anything to MT, I expect those to default
to a Monogenic MOI


DOES NOT CURRENTLY CHECK PARENT GENOTYPES - MVP holds that everyone is
a singleton
"""


import logging
from abc import abstractmethod
from typing import Any, Dict, List, Optional, Set

from reanalysis.pedigree import PedigreeParser
from reanalysis.utils import (
    AbstractVariant,
    CompHetDict,
    ReportedVariant,
)


# config keys to use for dominant MOI tests
GNOMAD_RARE_THRESHOLD = 'gnomad_dominant'
GNOMAD_AD_AC_THRESHOLD = 'gnomad_max_ac_dominant'
GNOMAD_DOM_HOM_THRESHOLD = 'gnomad_max_homs_dominant'
GNOMAD_REC_HOM_THRESHOLD = 'gnomad_max_homs_recessive'
INFO_HOMS = {'gnomad_hom', 'gnomad_ex_hom', 'exac_ac_hom'}


def check_for_second_hit(
    first_variant: str, comp_hets: CompHetDict, sample: str, gene: str
) -> List[str]:
    """
    checks for a second hit partner in this gene

    DOES NOT CURRENTLY CHECK VARIANT PHASE
    DOES NOT CURRENTLY CHECK PARENT GENOTYPES

    :param first_variant: string representation of variant1
    :param comp_hets: lookup dict for compound hets
    :param sample: ID string
    :param gene: gene string
    :return:
    """

    response = []

    # check if the sample has any comp-hets
    if sample not in comp_hets.keys():
        return response
    sample_dict = comp_hets.get(sample)

    # check if this sample-gene combination has listed variants
    if gene in sample_dict:

        # check if this variant has any listed partners in this gene
        gene_dict = sample_dict.get(gene)
        if first_variant in gene_dict:
            response = gene_dict.get(first_variant)

    return response


class MOIRunner:
    """
    pass
    """

    def __init__(
        self,
        pedigree: PedigreeParser,
        target_moi: str,
        config: Dict[str, Any],
        comp_het_lookup: CompHetDict,
    ):
        """
        for each possible MOI, choose the appropriate filters to apply
        ran into a situation where the ID of target_moi didn't match the
        exact same MOI as the IDs were different.

        This logic is only called once per MOI, not once per variant

        :param pedigree:
        :param target_moi:
        :param config:
        :param comp_het_lookup: the dictionary mapping all CH variant pairs
        """

        # for unknown, we catch all possible options?
        # should we be doing both checks for Monoallelic?
        if target_moi == 'Monoallelic':
            self.filter_list = [
                DominantAutosomal(pedigree=pedigree, config=config),
            ]
        elif target_moi in ['Mono_And_Biallelic', 'Unknown']:
            self.filter_list = [
                DominantAutosomal(pedigree=pedigree, config=config),
                RecessiveAutosomal(
                    pedigree=pedigree, config=config, comp_het=comp_het_lookup
                ),
            ]
        elif target_moi == 'Biallelic':
            self.filter_list = [
                RecessiveAutosomal(
                    pedigree=pedigree, config=config, comp_het=comp_het_lookup
                )
            ]

        elif target_moi == 'Hemi_Mono_In_Female':
            self.filter_list = [
                XRecessive(pedigree=pedigree, config=config, comp_het=comp_het_lookup),
                XDominant(pedigree=pedigree, config=config),
            ]

        elif target_moi == 'Hemi_Bi_In_Female':
            self.filter_list = [
                XRecessive(pedigree=pedigree, config=config, comp_het=comp_het_lookup)
            ]

        elif target_moi == 'Y_Chrom_Variant':
            self.filter_list = [YHemi(pedigree=pedigree, config=config)]

        else:
            raise Exception(f'MOI type {target_moi} is not addressed in MOI')

    def run(
        self, principal_var, gene_lookup: Dict[str, AbstractVariant]
    ) -> List[ReportedVariant]:
        """
        run method - triggers each relevant inheritance model
        :param principal_var: the variant we are focused on
        :param gene_lookup: all variants in the gene, indexed on chr-pos-ref-alt
        :return:
        """
        moi_matched = []
        for model in self.filter_list:
            moi_matched.extend(
                model.run(principal_var=principal_var, gene_lookup=gene_lookup)
            )
        return moi_matched

    def send_it(self):
        """
        to stop pylint complaining
        :return:
        """
        print(f'Yaaas {self.filter_list}')


class BaseMoi:
    """
    Definition of the MOI base class
    """

    def __init__(
        self,
        pedigree: PedigreeParser,
        config: Dict[str, Any],
        applied_moi: str,
        comp_het: Optional[CompHetDict],
    ):
        """
        base class
        """
        if applied_moi is None:
            raise Exception('An applied MOI needs to reach the Base Class')
        self.pedigree = pedigree
        self.config = config
        self.applied_moi = applied_moi
        self.comp_het = comp_het

    @abstractmethod
    def run(
        self, principal_var: AbstractVariant, gene_lookup: Dict[str, AbstractVariant]
    ) -> List[ReportedVariant]:
        """
        run all applicable inheritance patterns and finds good fits
        :param principal_var:
        :param gene_lookup:
        :return:
        """

    def is_affected(self, sample_id: str) -> bool:
        """
        take a sample ID and check if they are affected
        :param sample_id:
        :return:
        """
        if self.pedigree.participants[sample_id].details.affected:
            return True
        return False

    def check_familial_inheritance(
        self,
        sample_id: str,
        called_variants: Set[str],
        complete_penetrance: bool = True,
    ) -> bool:
        """
        sex-agnostic check for single variant inheritance
        designed to be called recursively, this method will take a single sample
        as a starting point and check that the MOI is viable across the whole family
        *we are assuming complete penetrance with this version*

        - find the family ID from this sample
        - iterate through all family members, and check that MOI holds for all

        Return False if a participant fails tests, True if all checked (so far) are ok

        :param sample_id:
        :param called_variants: the set of sample_ids which have this variant
        :param complete_penetrance: if True, (force affected==has variant call)

        NOTE: this called_variants pool is prepared before calling this method.
        If we only want to check hom calls, only send hom calls. If we are checking for
        dominant conditions, we would bundle hom and het calls into this set

        This check should work for both autosomal & sex chromosomes...
        for X-chrom biallelic, send female Homs, and all male variants
        For autosomal biallelic send male and female Homs
        :return:
        """

        # get the family ID for this test
        family_id = self.pedigree.participants[sample_id].details.family

        # initial value
        variant_passing: bool = True

        # iterate through all family members, no interested in directionality
        # of relationships at the moment
        for member in self.pedigree.families[family_id]:

            # complete & incomplete penetrance - affected samples must have the variant
            # complete pen. requires participants to be affected if they have the var
            # if any of these combinations occur, fail the family
            if (
                member.details.affected
                and member.details.sample_id not in called_variants
            ) or (
                member.details.sample_id in called_variants
                and complete_penetrance
                and not member.details.affected
            ):
                # fail
                return False

        return variant_passing

    def check_familial_comp_het(
        self,
        sample_id: str,
        called_variants_1: Set[str],
        called_variants_2: Set[str],
        complete_penetrance: bool = True,
    ) -> bool:
        """
        compound_het check, requires 2 pools of variant calls
        called_variants_1 & called_variants_2 are the relevant variant calls from
        2 different variants. We evaluate the sample as 'having this variant pair'
        if the sample ID appears in both call groups

        - find the family ID from this sample
        - iterate through all family members, and check that MOI holds for all
        - this MOI test requires the participant to have both variant calls in order
          to 'have this variant'

        At each stage, append the sample ID of all checked samples so we don't repeat
        One broken check will fail the variant for the whole family

        Return False if a participant fails tests, True if all members pass

        :param sample_id:
        :param called_variants_1: called samples for variant 1
        :param called_variants_2: called samples for variant 2
        :param complete_penetrance: if True, (force affected==has variant call)
        :return:
        """

        # get the family ID for this test
        family_id = self.pedigree.participants[sample_id].details.family

        # initial value
        variant_passing: bool = True

        # iterate through all family members, no interested in directionality
        # of relationships at the moment
        for member in self.pedigree.families[family_id]:

            # one value to store the check that this sample has _this_ comp-het
            sample_comp_het = (
                member.details.sample_id in called_variants_1
                and member.details.sample_id in called_variants_2
            )

            # complete & incomplete penetrance - affected samples must have the variant
            # complete pen. requires participants to be affected if they have the var
            # if any of these combinations occur, fail the family
            if (member.details.affected and not sample_comp_het) or (
                sample_comp_het and complete_penetrance and not member.details.affected
            ):
                # fail
                return False

        return variant_passing


class DominantAutosomal(BaseMoi):
    """
    This class can also be called by the X-linked Dominant, in which case the
    Applied_MOI by name is overridden
    """

    def __init__(
        self,
        pedigree: PedigreeParser,
        config: Dict[str, Any],
        applied_moi: str = 'Autosomal Dominant',
    ):
        """

        :param pedigree: not yet implemented
        :param config:
        :param applied_moi:
        """
        self.ad_threshold = config.get(GNOMAD_RARE_THRESHOLD)
        self.ac_threshold = config.get(GNOMAD_AD_AC_THRESHOLD)
        self.hom_threshold = config.get(GNOMAD_DOM_HOM_THRESHOLD)
        super().__init__(
            pedigree=pedigree, config=config, applied_moi=applied_moi, comp_het=None
        )

    def run(
        self, principal_var: AbstractVariant, gene_lookup: Dict[str, AbstractVariant]
    ) -> List[ReportedVariant]:
        """
        simplest
        if variant is present and sufficiently rare, we take it

        :param principal_var:
        :param gene_lookup:
        :return:
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
        ):
            return classifications

        # autosomal dominant doesn't require support, but consider het and hom
        samples_with_this_variant = principal_var.het_samples.union(
            principal_var.hom_samples
        )
        for sample_id in [
            sam for sam in samples_with_this_variant if self.is_affected(sam)
        ]:

            # check if this is a candidate for dominant inheritance
            if not self.check_familial_inheritance(
                sample_id=sample_id, called_variants=samples_with_this_variant
            ):
                continue

            classifications.append(
                ReportedVariant(
                    sample=sample_id,
                    gene=principal_var.info.get('gene_id'),
                    var_data=principal_var,
                    reasons={self.applied_moi},
                    supported=False,
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
        pedigree: PedigreeParser,
        config: Dict[str, Any],
        comp_het: CompHetDict,
        applied_moi: str = 'Autosomal Recessive',
    ):
        """ """
        self.hom_threshold = config.get(GNOMAD_REC_HOM_THRESHOLD)
        super().__init__(
            pedigree=pedigree, config=config, applied_moi=applied_moi, comp_het=comp_het
        )

    def run(
        self, principal_var: AbstractVariant, gene_lookup: Dict[str, AbstractVariant]
    ) -> List[ReportedVariant]:
        """
        valid if present as hom, or compound het
        counts as being phased if a compound het is split between parents
        Clarify if we want to consider a homozygous variant as 2 hets
        :param principal_var:
        :param gene_lookup:
        :return:
        """

        classifications = []

        # remove if too many homs are present in population databases
        # no stricter AF here - if we choose to, we can apply while labelling
        if any(
            {
                principal_var.info.get(hom_key, 0) >= self.hom_threshold
                for hom_key in INFO_HOMS
            }
        ):
            return classifications

        # homozygous is relevant directly
        for sample_id in [
            sam for sam in principal_var.hom_samples if self.is_affected(sam)
        ]:

            # check if this is a possible candidate for homozygous inheritance
            if not self.check_familial_inheritance(
                sample_id=sample_id, called_variants=principal_var.hom_samples
            ):
                continue

            classifications.append(
                ReportedVariant(
                    sample=sample_id,
                    gene=principal_var.info.get('gene_id'),
                    var_data=principal_var,
                    reasons={f'{self.applied_moi} Homozygous'},
                    supported=False,
                )
            )

        # if hets are present, try and find support
        for sample_id in [
            sam for sam in principal_var.het_samples if self.is_affected(sam)
        ]:

            for partner in check_for_second_hit(
                first_variant=principal_var.coords.string_format,
                comp_hets=self.comp_het,
                sample=sample_id,
                gene=principal_var.info.get('gene_id'),
            ):

                # check if this is a candidate for comp-het inheritance
                if not self.check_familial_comp_het(
                    sample_id=sample_id,
                    called_variants_1=principal_var.het_samples,
                    called_variants_2=gene_lookup[partner].het_samples,
                ):
                    continue

                classifications.append(
                    ReportedVariant(
                        sample=sample_id,
                        gene=principal_var.info.get('gene_id'),
                        var_data=principal_var,
                        reasons={f'{self.applied_moi} Compound-Het'},
                        supported=True,
                        support_vars=[partner],
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

    def __init__(
        self,
        pedigree: PedigreeParser,
        config: Dict[str, Any],
        applied_moi: str = 'X_Dominant',
    ):
        """
        accept male hets and homs, and female hets without support
        :param pedigree:
        :param config:
        :param applied_moi:
        """
        self.ad_threshold = config.get(GNOMAD_RARE_THRESHOLD)
        self.ac_threshold = config.get(GNOMAD_AD_AC_THRESHOLD)
        self.hom_threshold = config.get(GNOMAD_DOM_HOM_THRESHOLD)
        super().__init__(
            pedigree=pedigree, config=config, applied_moi=applied_moi, comp_het=None
        )

    def run(
        self, principal_var: AbstractVariant, gene_lookup: Dict[str, AbstractVariant]
    ) -> List[ReportedVariant]:
        """
        if variant is present and sufficiently rare, we take it

        :param principal_var:
        :param gene_lookup:
        :return:
        """
        classifications = []

        if principal_var.coords.chrom.lower() != 'x':
            raise Exception(
                f'X-Chromosome MOI given for variant on {principal_var.coords.chrom}'
            )

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
        ):
            return classifications

        # all samples which have a variant call
        samples_with_this_variant = principal_var.het_samples.union(
            principal_var.hom_samples
        )

        for sample_id in samples_with_this_variant:

            # check if this is a candidate for dominant inheritance
            if not self.check_familial_inheritance(
                sample_id=sample_id, called_variants=samples_with_this_variant
            ):
                continue

            # passed inheritance test, create the record
            sex = (
                'Female'
                if self.pedigree.participants[sample_id].details.is_female
                else 'Male'
            )
            classifications.append(
                ReportedVariant(
                    sample=sample_id,
                    gene=principal_var.info.get('gene_id'),
                    var_data=principal_var,
                    reasons={f'{self.applied_moi} {sex}'},
                    supported=False,
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
        pedigree: PedigreeParser,
        config: Dict[str, Any],
        comp_het: CompHetDict,
        applied_moi: str = 'X_Recessive',
    ):
        """
        set parameters specific to recessive tests
        and take the comp-het dictionary
        :param pedigree:
        :param config:
        :param comp_het:
        :param applied_moi:
        """

        self.hom_dom_threshold = config.get(GNOMAD_DOM_HOM_THRESHOLD)
        self.hom_rec_threshold = config.get(GNOMAD_REC_HOM_THRESHOLD)

        super().__init__(
            pedigree=pedigree, config=config, applied_moi=applied_moi, comp_het=comp_het
        )

    def run(
        self, principal_var: AbstractVariant, gene_lookup: Dict[str, AbstractVariant]
    ) -> List[ReportedVariant]:
        """

        :param principal_var:
        :param gene_lookup:
        :return:
        """

        if principal_var.coords.chrom.lower() != 'x':
            raise Exception(
                f'X-Chromosome MOI given for variant on {principal_var.coords.chrom}'
            )

        classifications = []

        # remove from analysis if too many homs are present in population databases
        if any(
            {
                principal_var.info.get(hom_key, 0) >= self.hom_dom_threshold
                for hom_key in INFO_HOMS
            }
        ):
            return classifications

        # X-relevant, we separate out male and females
        # combine het and hom here, we don't trust the variant callers
        males = {
            sam
            for sam in principal_var.het_samples.union(principal_var.hom_samples)
            if not self.pedigree.participants[sam].details.is_female
        }

        # split female calls into 2 categories
        het_females = {
            sam
            for sam in principal_var.het_samples
            if self.pedigree.participants[sam].details.is_female
        }
        hom_females = {
            sam
            for sam in principal_var.hom_samples
            if self.pedigree.participants[sam].details.is_female
        }

        # if het females are present, try and find support
        for sample_id in het_females:

            for partner in check_for_second_hit(
                first_variant=principal_var.coords.string_format,
                comp_hets=self.comp_het,
                sample=sample_id,
                gene=principal_var.info.get('gene_id'),
            ):

                # check if this is a candidate for comp-het inheritance
                # get all female het calls on the paired variant
                het_females_partner = {
                    sam
                    for sam in gene_lookup[partner].het_samples
                    if self.pedigree.participants[sam].details.is_female
                }
                if not self.check_familial_comp_het(
                    sample_id=sample_id,
                    called_variants_1=principal_var.het_samples,
                    called_variants_2=het_females_partner,
                ):
                    continue

                classifications.append(
                    ReportedVariant(
                        sample=sample_id,
                        gene=principal_var.info.get('gene_id'),
                        var_data=principal_var,
                        reasons={f'{self.applied_moi} Compound-Het Female'},
                        supported=True,
                        support_vars=[partner],
                    )
                )

        # remove from analysis if too many homs are present in population databases
        if any(
            {
                principal_var.info.get(hom_key, 0) >= self.hom_rec_threshold
                for hom_key in INFO_HOMS
            }
        ):
            return classifications

        # find all het males and hom females
        # assumption that the sample can only be hom if female?
        samples_to_check = males.union(hom_females)
        for sample_id in samples_to_check:

            # check if this is a possible candidate for homozygous inheritance
            if not self.check_familial_inheritance(
                sample_id=sample_id, called_variants=samples_to_check
            ):
                continue

            if not self.pedigree.participants[sample_id].details.is_female:
                reason = f'{self.applied_moi} Male'
            else:
                reason = f'{self.applied_moi} Female'
            classifications.append(
                ReportedVariant(
                    sample=sample_id,
                    gene=principal_var.info.get('gene_id'),
                    var_data=principal_var,
                    reasons={reason},
                    supported=False,
                )
            )
        return classifications


class YHemi(BaseMoi):
    """
    so... we should flag any female calls?
    otherwise treat as AD
    not expecting to use this
    """

    def __name__(self):
        return self.__name__()

    def __init__(
        self,
        pedigree: PedigreeParser,
        config: Dict[str, Any],
        applied_moi: str = 'Y_Hemi',
    ):
        """

        :param pedigree:
        :param config:
        :param applied_moi:
        """

        self.ad_threshold = config.get(GNOMAD_RARE_THRESHOLD)
        self.ac_threshold = config.get(GNOMAD_AD_AC_THRESHOLD)
        super().__init__(
            pedigree=pedigree, config=config, applied_moi=applied_moi, comp_het=None
        )

    def run(
        self, principal_var: AbstractVariant, gene_lookup: Dict[str, AbstractVariant]
    ) -> List[ReportedVariant]:
        """
        flag calls on Y which are Hom (maybe ok?) or female (bit weird)
        :param principal_var:
        :param gene_lookup:
        :return:
        """
        classifications = []

        # more stringent Pop.Freq checks for dominant
        if (
            principal_var.info.get('gnomad_af') >= self.ad_threshold
            or principal_var.info.get('gnomad_ac') >= self.ac_threshold
        ):
            return classifications

        # y chrom... called as hom? Shouldn't be possible
        # half expecting this with GATK...
        for sample_id in principal_var.hom_samples:
            logging.warning(f'Sample {sample_id} is a hom call on Y')

        # we don't expect any confident Y calls in females
        for sample_id in principal_var.het_samples.union(principal_var.hom_samples):
            if self.pedigree.participants[sample_id].details.is_female:
                logging.error(f'Sample {sample_id} is a female with call on Y')
            classifications.append(
                ReportedVariant(
                    sample=sample_id,
                    gene=principal_var.info.get('gene_id'),
                    var_data=principal_var,
                    reasons={self.applied_moi},
                    supported=False,
                )
            )

        return classifications


# create a DeNovo category?
