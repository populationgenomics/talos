#!/usr/bin/env python3


"""
runs between classification and publishing results
takes a number of inputs:
    - Classified VCF
    - JSON describing the present compound-het pairs
    - PanelApp data

reads in all compound het pairs
reads in all panelapp details
for each variant in each participant, check MOI in affected
participants relative to the MOI described in PanelApp
"""

import json
import logging
import sys
from collections import defaultdict
from datetime import datetime
from typing import Union

import click
from cyvcf2 import VCFReader
from peddy.peddy import Ped

from cpg_utils import to_path
from cpg_utils.config import get_config

from reanalysis.moi_tests import MOIRunner, PEDDY_AFFECTED
from reanalysis.utils import (
    canonical_contigs_from_vcf,
    get_cohort_config,
    filter_results,
    find_comp_hets,
    gather_gene_dict_from_contig,
    get_new_gene_map,
    read_json_from_path,
    CustomEncoder,
    GeneDict,
    ReportedVariant,
)

AMBIGUOUS_FLAG = 'Ambiguous Cat.1 MOI'
MALE_FEMALE = {'male', 'female'}


def set_up_moi_filters(
    panelapp_data: dict,
    pedigree: Ped,
) -> dict[str, MOIRunner]:
    """
    parse the panelapp data, and find all MOIs in this dataset
    for each unique MOI, set up a MOI filter instance
    save each one to a dictionary

    {MOI_string: MOI_runner (with a .run() method)}

    The MOI_runner uses a MOI string to select appropriate filters

    All logic regarding how MOI is applied, and which MOIs to
    apply to which PanelApp MOI descriptions is partitioned off into
    the MOI classes. All we need here is a Run() method, that returns
    either a list of results, or an empty list

    for every variant, we can then do a simple lookup using this
    dictionary to find the correct MOI runner, and run it
    that will return all matching MOIs for the variant

    This dictionary format means we only have to set up each once
    A billion variants, 6 MOI = 6 test instances, each created once

    Args:
        panelapp_data (dict):
        pedigree (Ped):

    Returns:
        a list
    """

    moi_dictionary = {}

    # iterate over all genes
    for gene_data in panelapp_data['genes'].values():

        # extract the per-gene MOI, don't re-simplify
        gene_moi = gene_data.get('moi')

        # if we haven't seen this MOI before, set up the appropriate filter
        if gene_moi not in moi_dictionary:
            # get a MOIRunner with the relevant filters
            moi_dictionary[gene_moi] = MOIRunner(pedigree=pedigree, target_moi=gene_moi)

    return moi_dictionary


def apply_moi_to_variants(
    variant_dict: GeneDict,
    moi_lookup: dict[str, MOIRunner],
    panelapp_data: dict[str, dict[str, Union[str, bool]]],
    pedigree: Ped,
) -> list[ReportedVariant]:
    """
    take all variants on a given contig & MOI filters
    find all variants/compound hets which fit the PanelApp MOI

    Args:
        variant_dict (dict): all possible variants, lists indexed by gene
        moi_lookup (dict): the MOI model runner per MOI string
        panelapp_data (dict): all genes and relevant details
        pedigree (Ped): the pedigree for this cohort
    """

    results = []

    for gene, variants in variant_dict.items():

        comp_het_dict = find_comp_hets(var_list=variants, pedigree=pedigree)

        # extract the panel data specific to this gene
        # extract once per gene, not once per variant
        panel_gene_data = panelapp_data.get(gene)

        # variant appears to be in a red gene
        if panel_gene_data is None:
            logging.error(f'How did this gene creep in? {gene}')
            continue

        for variant in variants:

            # is this even possible?
            if not (variant.het_samples or variant.hom_samples):
                continue

            # this variant is a candidate for MOI checks
            # - use MOI to get appropriate model
            # - run variant, append relevant classification(s) to the results
            # - always run partially penetrant analysis for Category 1 (clinvar)
            # pass on whether this variant is support only
            # - no dominant MOI
            # - discarded if two support-only form a comp-het
            variant_results = moi_lookup[panel_gene_data.get('moi')].run(
                principal_var=variant,
                comp_het=comp_het_dict,
                partial_pen=variant.info.get('categoryboolean1', False),
            )

            # Flag! If this is a Category 1 (ClinVar) variant, and we are
            # interpreting under a lenient MOI, add flag for analysts
            # control this in just one place
            if panel_gene_data.get('moi') == 'Mono_And_Biallelic' and variant.info.get(
                'categoryboolean1', False
            ):

                # consider each variant in turn
                for each_result in variant_results:

                    # never tag if this variant/sample is de novo
                    if '4' in each_result.var_data.categories:
                        continue

                    if each_result.reasons == {'Autosomal Dominant'}:
                        each_result.flags += [AMBIGUOUS_FLAG]

            results.extend(variant_results)

    return results


def clean_and_filter(
    results_holder: dict,
    result_list: list[ReportedVariant],
    panelapp_data: dict,
    participant_panels: dict | None = None,
) -> dict[str, list[ReportedVariant]]:
    """
    It's possible 1 variant can be classified multiple ways
    e.g. different MOIs (dominant and comp het)
    e.g. the same variant with annotation from two genes

    This cleans those to unique for final report
    stores panel names within the 'panels' attribute, either
    as matched (phenotype matched)
    as forced (cohort-wide applied panel)
    or neither

    Args:
        results_holder (): container for all results data
        result_list (): list of all ReportedVariant events
        panelapp_data ():
        participant_panels ():

    Returns:
        cleaned data
    """
    # pylint: disable=too-many-branches

    cohort_panels = set(get_cohort_config().get('cohort_panels', []))

    panel_meta = {
        content['id']: content['name'] for content in panelapp_data['metadata']
    }

    # if we have a lookup, grab the relevant information
    if participant_panels is not None:
        participant_panels = {
            sample: set(content['panels'])
            for sample, content in participant_panels.items()
        }

    gene_details = {}

    for each_event in result_list:

        each_event.independent = each_event.is_independent

        # grab some attributes from the event
        sample = each_event.sample
        gene = each_event.gene
        variant = each_event.var_data

        # no classifications = not interesting. Shouldn't be possible
        if not variant.categories:
            continue

        # find all panels for this gene
        if gene in gene_details:
            all_panels = gene_details[gene]

        else:
            # don't re-cast sets for every single variant
            all_panels = set(panelapp_data['genes'][gene]['panels'])
            gene_details[gene] = all_panels

        # get all forced panels this gene intersects with
        cohort_intersection: set = cohort_panels.intersection(all_panels)

        # check that the gene is in a panel of interest, and confirm new
        # neither step is required if no custom panel data is supplied
        if participant_panels is not None:

            # intersection to find participant phenotype-matched panels
            phenotype_intersection: set = participant_panels[sample].intersection(
                all_panels
            )

            # re-intersect to join phenotype matched with cohort-forced
            full_intersection = phenotype_intersection.union(cohort_intersection)

            # is this gene relevant for this participant?
            # this test includes matched, cohort-level, and core panel
            if not full_intersection:
                continue

            each_event.panels['matched'] = [
                panel_meta[pid]
                for pid in phenotype_intersection
                if pid != get_config()['workflow'].get('default_panel', 137)
            ]

        if cohort_intersection:
            each_event.panels['forced'] = [
                panel_meta[pid] for pid in cohort_intersection
            ]

        # equivalence logic might need a small change here -
        # If this variant and that variant have same sample/pos, equivalent
        # If either was independent, set that flag to True
        # Add a union of all Support Variants from both events
        if each_event not in results_holder[sample]['variants']:
            results_holder[sample]['variants'].append(each_event)

        else:
            prev_event = results_holder[sample]['variants'][
                results_holder[sample]['variants'].index(each_event)
            ]

            # if this is independent, set independent to True
            if each_event.independent:
                prev_event.independent = True

            # take the union of all supporting variants for both
            prev_event.support_vars.update(each_event.support_vars)

            prev_event.reasons.update(each_event.reasons)
            prev_event.gene = ','.join(prev_event.gene.split(',') + [each_event.gene])

            # combine flags across variants, and remove Ambiguous marking
            # if it's no longer appropriate
            both_flags = sorted(set(prev_event.flags + each_event.flags))
            if (
                prev_event.reasons != {'Autosomal Dominant'}
                and AMBIGUOUS_FLAG in both_flags
            ):
                both_flags.remove(AMBIGUOUS_FLAG)
            prev_event.flags = both_flags

    # organise the variants by chromosomal location
    for sample in results_holder:
        results_holder[sample]['variants'].sort()

    return results_holder


def count_families(pedigree: Ped, samples: list[str]) -> dict:
    """
    add metadata to results
    parsed during generation of the report
    most of these inputs aren't used...

    affected, male, female, and family sizes all at the same level
    maybe re-think this output structure for the report

    Args:
        pedigree (Ped): the Peddy pedigree object for the family
        samples (list): all the samples explicitly in this VCF

    Returns:
        A breakdown of all the family structures within this VCF
    """

    # contains all sample IDs for the given families
    family_dict = defaultdict(set)

    # the final dict of counts to return
    stat_counter = defaultdict(int)

    # iterate over samples in the VCF
    for sample_id in samples:
        ped_sample = pedigree[sample_id]
        family_dict[ped_sample.family_id].add(sample_id)

        # direct count of # each sex and # affected
        if ped_sample.sex in MALE_FEMALE:
            stat_counter[ped_sample.sex] += 1
        else:
            stat_counter['unknown_sex'] += 1
        if ped_sample.affected == PEDDY_AFFECTED:
            stat_counter['affected'] += 1

    # now count family sizes and structures
    for family_id, family_samples in family_dict.items():

        # bool flag - if we found a 'trio' don't also
        # count as family size 3
        trio_bool = False

        ped_family = pedigree.families[family_id]

        for trio in ped_family.trios():

            # check for a trio with all samples present
            if all(each.sample_id in family_samples for each in trio):
                trio_bool = True
                # if the proband has a sibling, call this a quad
                if list(trio[0].full_siblings):
                    stat_counter['quads'] += 1
                # otherwise a trio
                else:
                    stat_counter['trios'] += 1
                break

        # if we counted as a trio/quad, don't re-count
        if trio_bool:
            continue

        stat_counter[str(len(family_samples))] += 1

    return dict(stat_counter)


def prepare_results_shell(
    vcf_samples: list[str], pedigree: Ped, panel_data: dict | None, panelapp: dict
) -> dict:
    """
    prepare an empty dictionary for the results, feat. participant metadata
    Args:
        vcf_samples (): samples in the VCF header
        pedigree (): the Peddy PED object
        panel_data (): dictionary of per-participant panels, or None
        panelapp (): dictionary of gene data
    Returns:
        a pre-populated dict with sample metadata filled in
    """

    if panel_data is None:
        panel_data = {}

    # create an empty dict for all the samples
    sample_dict = {}

    ext_conf_path = get_config().get('dataset_specific', {}).get('external_lookup')
    external_map = read_json_from_path(ext_conf_path, default={})
    panel_meta = {content['id']: content['name'] for content in panelapp['metadata']}

    for sample in [
        sam.sample_id
        for sam in pedigree.samples()
        if sam.affected == PEDDY_AFFECTED and sam.sample_id in vcf_samples
    ]:
        family = pedigree.families[pedigree[sample].family_id]

        family_members = {
            member.sample_id: {
                'sex': str(member.sex)
                if str(member.sex) in {'male', 'female'}
                else 'unknown',
                'affected': member.affected == PEDDY_AFFECTED,
                'ext_id': external_map.get(member.sample_id, member.sample_id),
            }
            for member in family
        }
        sample_dict[sample] = {
            'variants': [],
            'metadata': {
                'ext_id': external_map.get(sample, sample),
                'family_id': pedigree[sample].family_id,
                'members': family_members,
                'phenotypes': panel_data.get(sample, {}).get('hpo_terms', []),
                'panel_ids': panel_data.get(sample, {}).get('panels', []),
                'panel_names': [
                    panel_meta[panel_id]
                    for panel_id in panel_data.get(sample, {}).get('panels', [])
                ],
            },
        }

    return sample_dict


@click.command
@click.option('--labelled_vcf', help='Category-labelled VCF')
@click.option('--out_json', help='Prefix to write JSON results to')
@click.option('--panelapp', help='Path to JSON file of PanelApp data')
@click.option('--pedigree', help='Path to joint-call PED file')
@click.option(
    '--input_path', help='source data', default='Not supplied', show_default=True
)
@click.option('--participant_panels', help='panels per participant', default=None)
def main(
    labelled_vcf: str,
    out_json: str,
    panelapp: str,
    pedigree: str,
    input_path: str,
    participant_panels: str | None = None,
):
    """
    VCFs used here should be small
    These have been pre-filtered to retain only a small number of candidate variants
    holding all the variants in memory should not be a challenge, no matter how large
    the cohort; if the variant number is large, the classes should be refined
    We expect approximately linear scaling with participants in the joint call

    Args:
        labelled_vcf (str): VCF output from Hail Labelling stage
        out_json (str): location to write output file
        panelapp (str): location of PanelApp data JSON
        pedigree (str): location of PED file
        input_path (str): VCF/MT used as input
        participant_panels (str): json of panels per participant
    """

    out_json = to_path(out_json)

    # parse the pedigree from the file
    pedigree_digest = Ped(pedigree)

    # parse panelapp data from dict
    panelapp_data = read_json_from_path(panelapp)

    # set up the inheritance checks
    moi_lookup = set_up_moi_filters(
        panelapp_data=panelapp_data, pedigree=pedigree_digest
    )

    # open the VCF using a cyvcf2 reader
    vcf_opened = VCFReader(labelled_vcf)

    participant_panels = read_json_from_path(participant_panels)

    # create the new gene map
    new_gene_map = get_new_gene_map(panelapp_data, participant_panels)

    result_list = []

    # obtain a set of all contigs with variants
    for contig in canonical_contigs_from_vcf(vcf_opened):
        # assemble {gene: [var1, var2, ..]}
        contig_dict = gather_gene_dict_from_contig(
            contig=contig,
            variant_source=vcf_opened,
            new_gene_map=new_gene_map,
            singletons=bool('singleton' in pedigree),
        )

        result_list.extend(
            apply_moi_to_variants(
                variant_dict=contig_dict,
                moi_lookup=moi_lookup,
                panelapp_data=panelapp_data['genes'],
                pedigree=pedigree_digest,
            )
        )

    # create a shell to store results in
    results_shell = prepare_results_shell(
        vcf_samples=vcf_opened.samples,
        pedigree=pedigree_digest,
        panel_data=participant_panels,
        panelapp=panelapp_data,
    )

    # remove duplicate and invalid variants
    analysis_results = clean_and_filter(
        results_holder=results_shell,
        result_list=result_list,
        panelapp_data=panelapp_data,
        participant_panels=participant_panels,
    )

    # annotate previously seen results using cumulative data file(s)
    analysis_results = filter_results(
        analysis_results, singletons=bool('singleton' in pedigree)
    )

    # create the full final output file
    final_results = {
        'results': analysis_results,
        'metadata': {
            'input_file': input_path,
            'cohort': get_config()['workflow']['dataset'],
            'run_datetime': get_config()['workflow'].get(
                'fake_date', f'{datetime.now():%Y-%m-%d %H:%M}'
            ),
            'family_breakdown': count_families(
                pedigree_digest, samples=vcf_opened.samples
            ),
            'panels': panelapp_data['metadata'],
            'container': get_config()['workflow']['driver_image'],
            'categories': get_config()['categories'],
        },
    }

    # store results using the custom-encoder to transform sets & DataClasses
    with out_json.open('w') as fh:
        json.dump(final_results, fh, cls=CustomEncoder, indent=4)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )
    main()  # pylint: disable=no-value-for-parameter
