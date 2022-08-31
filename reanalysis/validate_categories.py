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
from argparse import ArgumentParser
from collections import defaultdict
from typing import Any, Dict, List, Union

from cloudpathlib import AnyPath
from cyvcf2 import VCFReader
from peddy.peddy import Ped

from reanalysis.moi_tests import MOIRunner, PEDDY_AFFECTED
from reanalysis.utils import (
    canonical_contigs_from_vcf,
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
    config: Dict[str, Any],
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
    :param config:
    :param pedigree:
    :return:
    """

    moi_dictionary = {}

    # iterate over all genes
    for key, gene_data in panelapp_data.items():

        # skip over the stored metadata
        if '_version' in key:
            continue

        # extract the per-gene MOI, and SIMPLIFY
        gene_moi = get_simple_moi(gene_data.get('moi'))

        # if we haven't seen this MOI before, set up the appropriate filter
        if gene_moi not in moi_dictionary:

            # get a MOIRunner with the relevant filters
            moi_dictionary[gene_moi] = MOIRunner(
                pedigree=pedigree, target_moi=gene_moi, config=config['moi_tests']
            )

    return moi_dictionary


def apply_moi_to_variants(
    variant_dict: GeneDict,
    moi_lookup: Dict[str, MOIRunner],
    panelapp_data: Dict[str, Dict[str, Union[str, bool]]],
    pedigree: Ped,
) -> List[ReportedVariant]:
    """
    take a collection of all variants on a given contig & MOI filters
    find all variants/compound hets which fit the PanelApp MOI

    :param variant_dict:
    :param moi_lookup:
    :param panelapp_data:
    :param pedigree:
    :return:
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
    Possibility 1 variant can be classified multiple ways
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


def update_result_meta(
    results: dict, config: dict, pedigree: Ped, panelapp: dict
) -> dict:
    """
    takes the 'cleaned' results, and adds in a metadata key
    the key is used to set the results in context, and will be parsed
    during generation of the report
    """
    family_counter = defaultdict(int)
    for family in pedigree.families:
        affected, sex, trios, quads = pedigree.families[family].summary()
        family_counter['affected'] += affected[True]
        family_counter['male'] += sex['male']
        family_counter['female'] += sex['female']
        if trios != 0:
            family_counter['trios'] += trios
        if quads != 0:
            family_counter['quads'] += quads
        family_counter[str(len(pedigree.families[family].samples))] += 1

    panels = panelapp['metadata']['additional_panels']
    panels.append(
        {
            'panel_name': panelapp['metadata']['panel_name'],
            'panel_version': panelapp['metadata']['panel_version'],
        }
    )

    results['metadata'] = {
        'run_datetime': config['latest_run'],
        'input_file': config['input_file'],
        'family_breakdown': dict(family_counter),
        'panels': panels,
    }

    return results


def main(
    labelled_vcf: str,
    config_path: Union[str, Dict[str, Any]],
    out_json: str,
    panelapp: str,
    pedigree: str,
):
    """
    VCFs used here should be small
    These have been pre-filtered to retain only a small number of candidate variants
    holding all the variants in memory should not be a challenge, no matter how large
    the cohort; if the variant number is large, the classes should be refined
    We expect approximately linear scaling with participants in the joint call
    :param labelled_vcf:
    :param config_path:
    :param out_json:
    :param panelapp:
    :param pedigree:
    """

    # check if this is a singleton pedigree
    singletons = 'singleton' in pedigree

    # parse the pedigree from the file
    pedigree_digest = Ped(pedigree)

    # parse panelapp data from dict
    panelapp_data = read_json_from_path(panelapp)

    # get the runtime configuration
    if isinstance(config_path, dict):
        config_dict = config_path
    elif isinstance(config_path, str):
        config_dict = read_json_from_path(config_path)
    else:
        raise Exception(
            f'What is the conf path then?? "{config_path}": {type(config_path)}'
        )

    # set up the inheritance checks
    moi_lookup = set_up_inheritance_filters(
        panelapp_data=panelapp_data, pedigree=pedigree_digest, config=config_dict
    )

    # open the VCF using a cyvcf2 reader
    vcf_opened = VCFReader(labelled_vcf)

    # permit a blacklist to exclude known artefacts
    variant_blacklist = None
    if 'variant_blacklist' in config_dict:
        variant_blacklist = read_json_from_path(config_dict['variant_blacklist'])

    results = []

    # obtain a set of all contigs with variants
    for contig in canonical_contigs_from_vcf(vcf_opened):

        # assemble {gene: [var1, var2, ..]}
        contig_dict = gather_gene_dict_from_contig(
            contig=contig,
            variant_source=vcf_opened,
            config=config_dict,
            panelapp_data=panelapp_data,
            singletons=singletons,
            blacklist=variant_blacklist,
        )

        results.extend(
            apply_moi_to_variants(
                variant_dict=contig_dict,
                moi_lookup=moi_lookup,
                panelapp_data=panelapp_data,
                pedigree=pedigree_digest,
            )
        )

    # remove duplicate variants
    cleaned_results = clean_initial_results(
        results, samples=vcf_opened.samples, pedigree=pedigree_digest
    )

    # add metadata into the results
    meta_results = update_result_meta(
        cleaned_results, config_dict, pedigree=pedigree_digest, panelapp=panelapp_data
    )

    # dump results using the custom-encoder to transform sets & DataClasses
    with AnyPath(out_json).open('w') as fh:
        json.dump(meta_results, fh, cls=CustomEncoder, indent=4)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    parser = ArgumentParser()
    parser.add_argument('--labelled_vcf', help='Category-labelled VCF')
    parser.add_argument('--config_path', help='path to the runtime JSON config')
    parser.add_argument('--pedigree', help='Path to joint-call PED file')
    parser.add_argument('--panelapp', help='Path to JSON file of PanelApp data')
    parser.add_argument('--out_json', help='Path to write JSON results to')
    args = parser.parse_args()
    main(
        labelled_vcf=args.labelled_vcf,
        config_path=args.config_path,
        out_json=args.out_json,
        panelapp=args.panelapp,
        pedigree=args.pedigree,
    )
