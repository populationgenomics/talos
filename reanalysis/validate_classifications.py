"""
runs between classification and publishing results
takes 2 VCFs: classes and compound hets
reads in all compound het pairs
reads in all panelapp details
for each variant in each participant, check MOI
"""
import json
import logging
from typing import Any, Dict, List, Union
from functools import cache
import click
from cloudpathlib import AnyPath
from cyvcf2 import VCFReader

from reanalysis.moi_tests import MOIRunner
from reanalysis.utils import (
    AnalysisVariant,
    canonical_contigs_from_vcf,
    CompHetDict,
    CustomEncoder,
    get_simple_moi,
    parse_ped_simple,
    PanelAppDict,
    PedPerson,
    read_json_dict_from_path,
    ReportedVariant,
)


def set_up_inheritance_filters(
    panelapp_data: Dict[str, Dict[str, Union[str, bool]]],
    config: Dict[str, Any],
    pedigree: Dict[str, PedPerson],
    comp_het_lookup: CompHetDict,
) -> Dict[str, MOIRunner]:
    """
    parse the panelapp data, and find all MOIs in this dataset
    for each unique MOI, set up a MOI filter instance
    save each one to a dictionary

    {MOI_string: MOI_runner (with a .run() method)}

    The MOI_runner class will use the provided MOI string to
    select which filters will be appropriate

    All logic regarding how MOI is applied, and which MOIs to
    apply to which PanelApp MOI descriptions is partitioned off into
    the MOI class. All we need here is a Run() method, that returns
    either a list of results, or an empty list

    for every variant, we can then do a simple lookup using this
    dictionary to find the correct MOI runner, and run it
    that will return all matching MOIs for the variant

    This dictionary format means we only have to set up each once
    A billion variants, 6 MOI = 6 test instances, each created once
    :param panelapp_data:
    :param config:
    :param pedigree:
    :param comp_het_lookup:
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
                pedigree=pedigree,
                target_moi=gene_moi,
                config=config['moi_tests'],
                comp_het_lookup=comp_het_lookup,
            )

    return moi_dictionary


@cache
def get_moi_from_panelapp(
    panelapp_data: PanelAppDict, gene_name: str
) -> Union[Dict[str, str], None]:
    """

    :param panelapp_data:
    :param gene_name:
    :return:
    """

    # extract the panel data specific to this gene
    panel_gene_data = panelapp_data.get(gene_name)
    return panel_gene_data


# pylint: disable=too-many-locals
def apply_moi_to_variants(
    classified_variant_source: str,
    moi_lookup: Dict[str, MOIRunner],
    panelapp_data: Dict[str, Dict[str, Union[str, bool]]],
    config: Dict[str, Any],
) -> List[ReportedVariant]:
    """
    take the variant source and list of established MOI filters
    find all variants/compound hets for all samples which fit within
    the PanelApp provided MOI

    :param classified_variant_source:
    :param moi_lookup:
    :param panelapp_data:
    :param config:
    :return:
    """

    results = []

    variant_source = VCFReader(classified_variant_source)
    vcf_samples = variant_source.samples

    # split here - process each contig separately
    # once the variants are parsed into plain dicts (pickle-able)
    # we could run the MOI tests in parallel
    for contig in canonical_contigs_from_vcf(variant_source):

        # a dict to allow lookup of variants on this whole chromosome
        contig_variants = {}

        # iterate over all variants on this contig and store by unique key
        # if contig has no variants, prints an error and returns []
        for each_variant in variant_source(contig):
            analysis_var = AnalysisVariant(
                each_variant, samples=vcf_samples, config=config
            )
            contig_variants[analysis_var.coords.string_format] = analysis_var

        # NOTE - if we want details from both sides of a compound het in the
        # MOI check or report details, this lookup gives us a way to access that data
        # For now just iterate over the individual variants
        # crawl through all the results
        for variant in contig_variants.values():

            # each variant assoc. a single gene
            gene = variant.info.get('gene_id')

            # extract the panel data specific to this gene, with cache
            panel_gene_data = get_moi_from_panelapp(panelapp_data, gene)

            # variant appears to be in a red gene
            if panel_gene_data is None:
                logging.error(f'How did this gene creep in? {gene}')
                continue

            # initial implementation: we only accept class2 when the gene is New
            # for now this is covered in Hail, so this confirmation is redundant
            # -
            # if the gene isn't 'new' in PanelApp, remove Class2 flag
            # in future expand to the 'if MOI has changed' logic
            if variant.class_2 and not panel_gene_data.get('new', False):
                variant.class_2 = False

            # we never use a C4-only variant as a principal variant
            # and we don't consider a variant with no assigned classes
            #   - no classes shouldn't have reached this point
            if variant.class_4_only or not variant.is_classified:
                continue

            results.extend(
                moi_lookup[get_simple_moi(panelapp_data[gene].get('moi'))].run(
                    principal_var=variant, ensg=gene
                )
            )

    return results


def clean_initial_results(
    result_list: List[ReportedVariant],
) -> Dict[str, Dict[str, ReportedVariant]]:
    """
    Possibility 1 variant can be classified multiple ways
    This cleans those to unique for final report
    Join all possible classes for the condensed variants
    :param result_list:
    """

    clean_results: Dict[str, Dict[str, ReportedVariant]] = {}

    for each_instance in result_list:
        support_id = (
            ','.join(sorted(each_instance.support_vars))
            if each_instance.support_vars is not None
            else 'Unsupported'
        )
        var_uid = (
            f'{each_instance.var_data.coords.string_format}__'
            f'{each_instance.gene}__'
            f'{support_id}'
        )

        # create a section for this sample if it doesn't exist
        if each_instance.sample not in clean_results:
            clean_results[each_instance.sample] = {}

        # get an existing object, or use the current one
        variant_object = clean_results.setdefault(each_instance.sample, {}).setdefault(
            var_uid, each_instance
        )

        # combine any possible reasons, and add
        clean_results[each_instance.sample][
            var_uid
        ].reasons = variant_object.reasons.union(each_instance.reasons)
    return clean_results


@click.command()
@click.option('--config_path', help='Path to a config JSON file')
@click.option(
    '--class_vcf',
    help='VCF from Hail with variant categories and VEP annotations',
)
@click.option('--comp_het', help='JSON containing all comp-het pairings')
@click.option('--pedigree', help='Pedigree file')
@click.option('--panelapp', help='PanelApp JSON file')
@click.option('--out_json', help='Write the analysis results in JSON form')
def main(
    config_path: str,
    class_vcf: str,
    comp_het: str,
    pedigree: str,
    panelapp: str,
    out_json: str,
):
    """
    All VCFs in use at this point will be small
    These have been pre-filtered to retain only a small number of classified variants
    holding all the variants in memory should not be a challenge, no matter how large
    the cohort; if the variant number is large, the classes should be refined
    We expect approximately linear scaling with participants in the joint call

    Might be able to use a single output path, just altering the extension
    Depends on how this is handled by Hail, as the object paths are Resource File paths

    Re-working of the comp-het logic means that we only store pairings as strings
    Not needing to reach the annotations attached to variant pairs opens up choices:
        - process each variant in turn (original design)
        - parse each chromosome separately, then process the group of variants
        - parse all variants, then process as a group

    these come with incrementing memory footprints...

    preference is for #2; process an entire contig together
        - we can look-up the partner variant's attributes in future if we want
        - this will be required when we do familial checks, e.g. for a compound het,
            we need the pair of variants to be absent in unaffected family members

    :param config_path:
    :param class_vcf:
    :param comp_het:
    :param pedigree:
    :param panelapp:
    :param out_json:
    """

    # parse the pedigree from the file
    pedigree_digest = parse_ped_simple(pedigree)

    # parse panelapp data from dict
    panelapp_data = read_json_dict_from_path(panelapp)

    # get the runtime configuration
    config_dict = read_json_dict_from_path(config_path)

    # find all the Compound Hets from C-H JSON
    comp_het_digest: CompHetDict = read_json_dict_from_path(comp_het)

    # set up the inheritance checks
    moi_lookup = set_up_inheritance_filters(
        panelapp_data=panelapp_data,
        pedigree=pedigree_digest,
        config=config_dict,
        comp_het_lookup=comp_het_digest,
    )

    # find classification events
    results = apply_moi_to_variants(
        classified_variant_source=class_vcf,
        moi_lookup=moi_lookup,
        panelapp_data=panelapp_data,
        config=config_dict,
    )

    # remove duplicate variants
    cleaned_results = clean_initial_results(results)

    # dump the JSON-results to a file, serialised & with AnyPath
    # use the custom-encoder to print sets and DataClasses
    serialised = json.dumps(cleaned_results, cls=CustomEncoder, indent=4)
    AnyPath(out_json).write_text(serialised)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=E1120
