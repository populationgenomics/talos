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

from collections import defaultdict

import click
from cpg_utils import to_path
from cpg_utils.config import get_config
from cyvcf2 import VCFReader
from peddy.peddy import Ped

from reanalysis.models import (
    CATEGORY_DICT,
    FamilyMembers,
    PanelApp,
    PanelDetail,
    ParticipantHPOPanels,
    ParticipantMeta,
    ParticipantResults,
    PhenotypeMatchedPanels,
    ReportVariant,
    ResultData,
    ResultMeta,
    ReportPanel,
)
from reanalysis.moi_tests import MOIRunner, PEDDY_AFFECTED
from reanalysis.static_values import get_logger
from reanalysis.utils import (
    canonical_contigs_from_vcf,
    filter_results,
    find_comp_hets,
    gather_gene_dict_from_contig,
    get_cohort_config,
    get_new_gene_map,
    read_json_from_path,
    GeneDict,
)

AMBIGUOUS_FLAG = 'Ambiguous Cat.1 MOI'
MALE_FEMALE = {'male', 'female'}


def set_up_moi_filters(
    panelapp_data: PanelApp,
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
        panelapp_data (PanelApp):
        pedigree (Ped):

    Returns:
        a dict of all MOI classes, indexed by MOI string
    """

    moi_dictionary = {}

    # iterate over all genes
    for gene_data in panelapp_data.genes.values():

        # extract the per-gene MOI, don't re-simplify
        gene_moi = gene_data.moi

        # if we haven't seen this MOI before, set up the appropriate filter
        if gene_moi not in moi_dictionary:
            # get a MOIRunner with the relevant filters
            moi_dictionary[gene_moi] = MOIRunner(pedigree=pedigree, target_moi=gene_moi)

    return moi_dictionary


def apply_moi_to_variants(
    variant_dict: GeneDict,
    moi_lookup: dict[str, MOIRunner],
    panelapp_data: dict[str, PanelDetail],
    pedigree: Ped,
) -> list[ReportVariant]:
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
            get_logger().error(f'How did this gene creep in? {gene}')
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
            panel_moi = panel_gene_data.moi
            runner = moi_lookup[panel_moi]
            assert isinstance(runner, MOIRunner)
            variant_results = runner.run(
                principal_var=variant,
                comp_het=comp_het_dict,
                partial_pen=bool(variant.info.get('categoryboolean1', False)),
            )

            # Flag! If this is a Category 1 (ClinVar) variant, and we are
            # interpreting under a lenient MOI, add flag for analysts
            # control this in just one place
            if panel_gene_data.moi == 'Mono_And_Biallelic' and variant.info.get(
                'categoryboolean1', False
            ):

                # consider each variant in turn
                for each_result in variant_results:

                    # never tag if this variant/sample is de novo
                    if '4' in each_result.categories:
                        continue

                    if each_result.reasons == {'Autosomal Dominant'}:
                        each_result.flags.add(AMBIGUOUS_FLAG)

            results.extend(variant_results)

    return results


def clean_and_filter(
    results_holder: ResultData,
    result_list: list[ReportVariant],
    panelapp_data: PanelApp,
    dataset: str,
    participant_panels: PhenotypeMatchedPanels | None = None,
) -> ResultData:
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
        result_list (): list of all ReportVariant events
        panelapp_data ():
        dataset (str): dataset to use for getting the config portion
        participant_panels ():

    Returns:
        cleaned data
    """
    cohort_panels = set(get_cohort_config(dataset).get('cohort_panels', []))

    panel_meta: dict[int, str] = {
        content.id: content.name for content in panelapp_data.metadata
    }

    gene_details: dict[str, set[int]] = {}

    for each_event in result_list:

        # shouldn't be possible, here as a precaution
        assert (
            each_event.categories
        ), f'No categories for {each_event.var_data.coordinates.string_format}'

        # find all panels for this gene
        if each_event.gene in gene_details:
            all_panels = gene_details[each_event.gene]

        else:
            # don't re-cast sets for every single variant
            all_panels = set(panelapp_data.genes[each_event.gene].panels)
            gene_details[each_event.gene] = all_panels

        # get all forced panels this gene intersects with
        cohort_intersection: set[int] = cohort_panels.intersection(all_panels)

        matched_panels = set()
        # check that the gene is in a panel of interest, and confirm new
        # neither step is required if no custom panel data is supplied
        if participant_panels is not None:

            # intersection to find participant phenotype-matched panels
            phenotype_intersection = participant_panels.samples[
                each_event.sample
            ].panels.intersection(all_panels)

            # is this gene relevant for this participant?
            # this test includes matched, cohort-level, and core panel
            if not phenotype_intersection.union(cohort_intersection):
                continue

            matched_panels = {
                panel_meta[pid]
                for pid in phenotype_intersection
                if pid != get_config()['workflow'].get('default_panel', 137)
            }

        forced_panels = set()
        if cohort_intersection:
            forced_panels = {panel_meta[pid] for pid in cohort_intersection}

        each_event.panels = ReportPanel(matched=matched_panels, forced=forced_panels)

        # equivalence logic might need a small change here -
        # If this variant and that variant have same sample/pos, equivalent
        # If either was independent, set that flag to True
        # Add a union of all Support Variants from both events
        if each_event not in results_holder.results[each_event.sample].variants:
            results_holder.results[each_event.sample].variants.append(each_event)

        else:
            prev_event = results_holder.results[each_event.sample].variants[
                results_holder.results[each_event.sample].variants.index(each_event)
            ]

            # if this is independent, set independent to True
            if each_event.independent:
                prev_event.independent = True

            # take the union of all supporting variants for both
            prev_event.support_vars.update(each_event.support_vars)

            prev_event.reasons.update(each_event.reasons)
            prev_event.gene = ','.join({*prev_event.gene.split(','), each_event.gene})

            # combine flags across variants, and remove Ambiguous marking
            # if it's no longer appropriate
            both_flags = {*prev_event.flags, *each_event.flags}
            if (
                prev_event.reasons != {'Autosomal Dominant'}
                and AMBIGUOUS_FLAG in both_flags
            ):
                both_flags.remove(AMBIGUOUS_FLAG)
            prev_event.flags = both_flags

    # organise the variants by chromosomal location... why?
    for sample in results_holder.results:
        results_holder.results[sample].variants.sort()

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
    family_dict: dict[str, set[str]] = defaultdict(set)

    # the final dict of counts to return
    stat_counter: dict[str, int] = defaultdict(int)

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
    results_meta: ResultMeta,
    vcf_samples: list[str],
    pedigree: Ped,
    dataset: str,
    panelapp: PanelApp,
    panel_data: PhenotypeMatchedPanels | None = None,
) -> ResultData:
    """
    Creates a ResultData object, with participant metadata filled out

    Args:
        results_meta (): metadata for the results
        vcf_samples (): samples in the VCF header
        pedigree (): the Peddy PED object
        dataset (str): dataset to use for getting the config portion
        panel_data (): dictionary of per-participant panels, or None
        panelapp (): dictionary of gene data

    Returns:
        ResultData with sample metadata filled in
    """

    if panel_data is None:
        panel_data = PhenotypeMatchedPanels()

    # create an empty dict for all the samples
    results_shell = ResultData(metadata=results_meta)

    # find the solved cases in this project
    solved_cases = get_cohort_config(dataset).get('solved_cases', [])
    panel_meta = {content.id: content.name for content in panelapp.metadata}

    # limit to affected samples present in both Pedigree and VCF
    for sample in [
        sam
        for sam in pedigree.samples()
        if sam.affected == PEDDY_AFFECTED and sam.sample_id in vcf_samples
    ]:
        sample_id = sample.sample_id
        family_id = sample.family_id

        family_members = {
            member.sample_id: FamilyMembers(
                **{
                    'sex': str(member.sex)
                    if str(member.sex) in {'male', 'female'}
                    else 'unknown',
                    'affected': member.affected == PEDDY_AFFECTED,
                    'ext_id': panel_data.samples.get(
                        member.sample_id, ParticipantHPOPanels()
                    ).external_id
                    or member.sample_id,
                }
            )
            for member in pedigree.families[family_id]
        }
        sample_panel_data = panel_data.samples.get(sample_id, ParticipantHPOPanels())
        results_shell.results[sample_id] = ParticipantResults(
            **{
                'variants': [],
                'metadata': ParticipantMeta(
                    **{
                        'ext_id': sample_panel_data.external_id or sample_id,
                        'family_id': pedigree[sample_id].family_id,
                        'members': family_members,
                        'phenotypes': sample_panel_data.hpo_terms,
                        'panel_ids': sample_panel_data.panels,
                        'panel_names': [
                            panel_meta[panel_id]
                            for panel_id in sample_panel_data.panels
                        ],
                        'solved': bool(
                            sample_id in solved_cases or family_id in solved_cases
                        ),
                    }
                ),
            }
        )

    return results_shell


@click.command
@click.option('--labelled_vcf', help='Category-labelled VCF')
@click.option('--labelled_sv', help='Category-labelled SV VCF', default=None)
@click.option('--out_json', help='Prefix to write JSON results to')
@click.option('--panelapp', help='Path to JSON file of PanelApp data')
@click.option('--pedigree', help='Path to joint-call PED file')
@click.option(
    '--input_path', help='source data', default='Not supplied', show_default=True
)
@click.option('--participant_panels', help='panels per participant', default=None)
@click.option('--dataset', help='optional, dataset to use', default=None)
def main(
    labelled_vcf: str,
    out_json: str,
    panelapp: str,
    pedigree: str,
    input_path: str,
    labelled_sv: str | None = None,
    participant_panels: str | None = None,
    dataset: str | None = None,
):
    """
    VCFs used here should be small
    These have been pre-filtered to retain only a small number of candidate variants
    holding all the variants in memory should not be a challenge, no matter how large
    the cohort; if the variant number is large, the classes should be refined
    We expect approximately linear scaling with participants in the joint call

    Args:
        labelled_vcf (str): VCF output from Hail Labelling stage
        labelled_sv (str | None): optional second VCF (SV)
        out_json (str): location to write output file
        panelapp (str): location of PanelApp data JSON
        pedigree (str): location of PED file
        input_path (str): VCF/MT used as input
        participant_panels (str): json of panels per participant
        dataset (str): optional, dataset to use
    """

    out_json_path = to_path(out_json)

    # parse the pedigree from the file
    ped = Ped(pedigree)

    # parse panelapp data from dict
    panelapp_data: PanelApp = read_json_from_path(panelapp, return_model=PanelApp)  # type: ignore

    # set up the inheritance checks
    moi_lookup = set_up_moi_filters(panelapp_data=panelapp_data, pedigree=ped)

    pheno_panels: PhenotypeMatchedPanels | None = read_json_from_path(
        participant_panels, return_model=PhenotypeMatchedPanels, default=None  # type: ignore
    )

    # create the new gene map
    new_gene_map = get_new_gene_map(panelapp_data, pheno_panels, dataset)

    result_list: list[ReportVariant] = []

    # open the small variant VCF using a cyvcf2 reader
    vcf_opened = VCFReader(labelled_vcf)

    # optional SV behaviour
    sv_opened = VCFReader(labelled_sv) if labelled_sv else None

    # obtain a set of all contigs with variants
    for contig in canonical_contigs_from_vcf(vcf_opened):
        # assemble {gene: [var1, var2, ..]}
        contig_dict = gather_gene_dict_from_contig(
            contig=contig,
            variant_source=vcf_opened,
            sv_source=sv_opened,
            new_gene_map=new_gene_map,
            singletons=bool('singleton' in pedigree),
        )

        result_list.extend(
            apply_moi_to_variants(
                variant_dict=contig_dict,
                moi_lookup=moi_lookup,
                panelapp_data=panelapp_data.genes,
                pedigree=ped,
            )
        )

    # create the full final output file
    results_meta = ResultMeta(
        **{
            'input_file': input_path,
            'cohort': dataset or get_config()['workflow']['dataset'] or 'unknown',
            'family_breakdown': count_families(ped, samples=vcf_opened.samples),
            'panels': panelapp_data.metadata,
            'container': get_config()['workflow']['driver_image'],
        }
    )

    # create a shell to store results in, adds participant metadata
    results_model = prepare_results_shell(
        results_meta=results_meta,
        vcf_samples=vcf_opened.samples,
        pedigree=ped,
        panel_data=pheno_panels,
        dataset=dataset,
        panelapp=panelapp_data,
    )

    # remove duplicate and invalid variants
    results_model = clean_and_filter(
        results_holder=results_model,
        result_list=result_list,
        panelapp_data=panelapp_data,
        dataset=dataset,
        participant_panels=pheno_panels,
    )

    # annotate previously seen results using cumulative data file(s)
    filter_results(
        results_model, singletons=bool('singleton' in pedigree), dataset=dataset
    )

    # write the output to long term storage using Pydantic
    # validate the model against the schema, then write the result if successful
    with to_path(out_json_path).open('w') as out_file:
        out_file.write(
            ResultData.model_validate(results_model).model_dump_json(indent=4)
        )


if __name__ == '__main__':
    get_logger(__file__).info('Starting MOI testing phase')
    get_logger().info(f'Operational Categories: {CATEGORY_DICT}')
    main()
