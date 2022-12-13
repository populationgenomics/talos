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
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Union

import click

from cyvcf2 import VCFReader
from peddy.peddy import Ped

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.git import get_git_commit_ref_of_current_repository

from reanalysis.moi_tests import MOIRunner, PEDDY_AFFECTED
from reanalysis.utils import (
    canonical_contigs_from_vcf,
    filter_results,
    find_comp_hets,
    gather_gene_dict_from_contig,
    get_simple_moi,
    read_json_from_path,
    CustomEncoder,
    GeneDict,
    ReportedVariant,
)


def set_up_inheritance_filters(
    panelapp_data: Dict[str, Dict[str, Union[str, bool]]],
    pedigree: Ped,
) -> Dict[str, MOIRunner]:
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
    :param panelapp_data:
    :param pedigree:
    :return:
    """

    moi_dictionary = {}

    # iterate over all genes
    for key, gene_data in panelapp_data.items():

        if key == 'metadata':
            continue

        # extract the per-gene MOI, and SIMPLIFY
        gene_moi = get_simple_moi(gene_data.get('moi'))

        # if we haven't seen this MOI before, set up the appropriate filter
        if gene_moi not in moi_dictionary:

            # get a MOIRunner with the relevant filters
            moi_dictionary[gene_moi] = MOIRunner(pedigree=pedigree, target_moi=gene_moi)

    return moi_dictionary


def apply_moi_to_variants(
    variant_dict: GeneDict,
    moi_lookup: Dict[str, MOIRunner],
    panelapp_data: Dict[str, Dict[str, Union[str, bool]]],
    pedigree: Ped,
) -> List[ReportedVariant]:
    """
    take all variants on a given contig & MOI filters
    find all variants/compound hets which fit the PanelApp MOI

    Args:
        variant_dict ():
        moi_lookup ():
        panelapp_data ():
        pedigree ():

    Returns:

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

        simple_moi = get_simple_moi(panel_gene_data.get('moi'))
        additional_panels = panel_gene_data.get('flags', [])

        for variant in variants:

            if not (variant.het_samples or variant.hom_samples):
                continue

            # if this variant is category 1, 2, 3, or 4; evaluate is as a 'primary'
            if variant.category_non_support:

                # this variant is a candidate for MOI checks
                # - find the simplified MOI string
                # - use to get appropriate MOI model
                # - run variant, append relevant classification(s) to the results
                # NEW - run partially penetrant analysis for Category 1 (clinvar)
                # adds a flag extension to include any specific panels for this gene
                variant_reports = moi_lookup[simple_moi].run(
                    principal_var=variant,
                    comp_het=comp_het_dict,
                    partial_penetrance=variant.info.get('categoryboolean1', False),
                )
                for var in variant_reports:
                    var.flags.extend(additional_panels)

                results.extend(variant_reports)

    return results


def clean_initial_results(
    result_list: list[ReportedVariant], samples: list[str], pedigree: Ped
) -> dict[str, list[ReportedVariant]]:
    """
    It's possible 1 variant can be classified multiple ways
    This cleans those to unique for final report
    Join all possible classes for the condensed variants
    :param result_list:
    :param samples: all samples from the VCF
    :param pedigree:
    """

    clean_results = defaultdict(list)

    for each_event in result_list:
        if each_event not in clean_results[each_event.sample]:
            clean_results[each_event.sample].append(each_event)
        else:
            prev_event_index = clean_results[each_event.sample].index(each_event)
            prev_event = clean_results[each_event.sample][prev_event_index]
            prev_event.reasons.update(each_event.reasons)
            prev_genes = set(prev_event.gene.split(','))
            prev_genes.add(each_event.gene)
            prev_event.gene = ','.join(prev_genes)
            prev_event.flags.extend(each_event.flags)
            prev_event.flags = list(set(prev_event.flags))

    # organise the variants by chromosomal location
    for sample in clean_results:
        clean_results[sample].sort()

    # Empty list for 0 variant samples with affected status
    # explicitly record samples checked in this analysis
    # the PED could have more samples than a joint call, due to sub-setting or QC.
    # When presenting results, we want all samples with negative findings, without
    # comparing both VCF and PED files
    affected_samples = [
        sam.sample_id
        for sam in pedigree.samples()
        if sam.affected == PEDDY_AFFECTED and sam.sample_id in samples
    ]
    for sample in affected_samples:
        if sample not in clean_results:
            clean_results[sample] = []

    return clean_results


def get_gene_panel_sets(gene_details: dict, gene: str) -> tuple[set, set]:
    """
    get each gene's associated panels only once

    shove in some lru_cache'ing here, so we don't keep generating the sets

    Args:
        gene_details ():
        gene ():

    Returns:
        set of all panels for this gene,
        set of new panels for this gene
    """
    single_gene_details = gene_details['genes'][gene]
    all_panels = set(single_gene_details['panels'])
    new_panels = set(single_gene_details['new'])
    return all_panels, new_panels


def gene_clean_results(
    party_panels: dict,
    panel_app_data: dict,
    cleaned_data: dict[str, list[ReportedVariant]],
) -> dict:
    """
    takes the unique'd data from the previous cleaning
    applies gene panel filters per-participant
    the only remaining data should be relevant to the participant's panel list

    This might be a bit heavy, but very few variants remain so runtime is meh

    Args:
        party_panels (): JSON file containing the panels per participant
        panel_app_data (): all panelapp results
        cleaned_data (): reportable variants, indexed per participant

    Returns:
        a gene-list filtered version of the reportable data
    """

    gene_cleaned_data = defaultdict(list)

    # panel lookup index is CPG ID
    for sample, variants in cleaned_data.items():

        participant_panels = set(party_panels[sample]['panels'])
        for variant in variants:
            all_panels, new_panels = get_gene_panel_sets(panel_app_data, variant.gene)

            if not bool(participant_panels.intersection(all_panels)):
                # this gene is not relevant for this participant
                continue

            # now check if Cat 2 is present - requires new gene status
            # pop it out if cat2 doesn't apply for this participant
            if '2' in variant.var_data.categories and not bool(
                participant_panels.intersection(new_panels)
            ):
                _ = variant.var_data.categories.pop(
                    variant.var_data.categories.index('2')
                )
                # should not be treated as new
                logging.info(f'Removing category 2 in {variant.gene} for {sample}')

            # if categories still remain, add the variant
            if variant.var_data.categories:
                gene_cleaned_data[sample].append(variant)

    return dict(gene_cleaned_data)


def count_families(pedigree: Ped, samples: list[str]) -> dict:
    """
    add metadata to results
    parsed during generation of the report
    most of these inputs aren't used...

    affected, male, female, and family sizes all at the same level
    maybe re-think this output structure for the report

    Args:
        pedigree ():
        samples ():

    Returns:

    """
    family_counter = defaultdict(int)
    for family in pedigree.families:
        # don't count families who don't appear in this pedigree subset
        if not any(sam.sample_id in samples for sam in pedigree.families[family]):
            continue

        affected, sex, trios, quads = pedigree.families[family].summary()
        family_counter['affected'] += affected[True]
        family_counter['male'] += sex['male']
        family_counter['female'] += sex['female']
        if trios != 0:
            family_counter['trios'] += trios
        if quads != 0:
            family_counter['quads'] += quads
        family_counter[str(len(pedigree.families[family].samples))] += 1

    return dict(family_counter)


@click.command
@click.option('--labelled_vcf', help='Category-labelled VCF')
@click.option('--out_json', help='Prefix to write JSON results to')
@click.option('--panelapp', help='Path to JSON file of PanelApp data')
@click.option('--pedigree', help='Path to joint-call PED file')
@click.option(
    '--input_path', help='source data', default='Not supplied', show_default=True
)
@click.option(
    '--participant_panels',
    help='dict of panels per participant',
    default=None,
)
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
    :param labelled_vcf:
    :param out_json:
    :param panelapp:
    :param pedigree:
    :param input_path: data file used as input
    :param participant_panels: the json of panels per participant
    """

    # parse the pedigree from the file
    pedigree_digest = Ped(pedigree)

    # parse panelapp data from dict
    panelapp_data = read_json_from_path(panelapp)

    # set up the inheritance checks
    moi_lookup = set_up_inheritance_filters(
        panelapp_data=panelapp_data, pedigree=pedigree_digest
    )

    # open the VCF using a cyvcf2 reader
    vcf_opened = VCFReader(labelled_vcf)

    results = []

    # obtain a set of all contigs with variants
    for contig in canonical_contigs_from_vcf(vcf_opened):

        # assemble {gene: [var1, var2, ..]}
        contig_dict = gather_gene_dict_from_contig(
            contig=contig,
            variant_source=vcf_opened,
            panelapp_data=panelapp_data,
            singletons=bool('singleton' in pedigree),
        )

        results.extend(
            apply_moi_to_variants(
                variant_dict=contig_dict,
                moi_lookup=moi_lookup,
                panelapp_data=panelapp_data,
                pedigree=pedigree_digest,
            )
        )

    # remove duplicate variants (better solution pls)
    analysis_results = clean_initial_results(
        results, samples=vcf_opened.samples, pedigree=pedigree_digest
    )

    # remove previously seen results using cumulative data files
    analysis_results = filter_results(
        analysis_results, singletons=bool('singleton' in pedigree)
    )

    # do we need to do multi-panel filtering?
    if participant_panels:
        analysis_results = gene_clean_results(
            party_panels=read_json_from_path(participant_panels),
            panel_app_data=panelapp_data,
            cleaned_data=analysis_results,
        )

    final_results = {
        'results': analysis_results,
        'metadata': {
            'input_file': input_path,
            'cohort': get_config()['workflow']['dataset'],
            'run_datetime': f'{datetime.now():%Y-%m-%d %H:%M}',
            'family_breakdown': count_families(
                pedigree_digest, samples=vcf_opened.samples
            ),
            'panels': panelapp_data['metadata'],
            'commit_id': get_git_commit_ref_of_current_repository(),
        },
    }

    # store results using the custom-encoder to transform sets & DataClasses
    with to_path(out_json).open('w') as fh:
        json.dump(final_results, fh, cls=CustomEncoder, indent=4)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=no-value-for-parameter
