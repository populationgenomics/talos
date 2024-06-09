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

from argparse import ArgumentParser
from collections import defaultdict

from cyvcf2 import VCFReader

from cpg_utils import to_path
from cpg_utils.config import config_retrieve

from talos.models import (
    CATEGORY_DICT,
    FamilyMembers,
    PanelApp,
    PanelDetail,
    ParticipantHPOPanels,
    ParticipantMeta,
    ParticipantResults,
    Pedigree,
    PhenotypeMatchedPanels,
    ReportPanel,
    ReportVariant,
    ResultData,
    ResultMeta,
)
from talos.moi_tests import MOIRunner
from talos.static_values import get_logger
from talos.utils import (
    GeneDict,
    canonical_contigs_from_vcf,
    filter_results,
    find_comp_hets,
    gather_gene_dict_from_contig,
    get_cohort_config,
    get_cohort_seq_type_conf,
    get_new_gene_map,
    make_flexible_pedigree,
    read_json_from_path,
)

AMBIGUOUS_FLAG = 'Ambiguous Cat.1 MOI'
MALE_FEMALE = {'1': 'male', '2': 'female', '-9': 'unknown_sex', '0': 'unknown_sex'}


def set_up_moi_filters(panelapp_data: PanelApp, pedigree: Pedigree) -> dict[str, MOIRunner]:
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
        pedigree (Pedigree):

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
    pedigree: Pedigree,
) -> list[ReportVariant]:
    """
    take all variants on a given contig & MOI filters
    find all variants/compound hets which fit the PanelApp MOI

    Args:
        variant_dict (dict): all possible variants, lists indexed by gene
        moi_lookup (dict): the MOI model runner per MOI string
        panelapp_data (dict): all genes and relevant details
        pedigree (Pedigree): the pedigree for this cohort
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
            if panel_gene_data.moi == 'Mono_And_Biallelic' and variant.info.get('categoryboolean1', False):
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

    - addition! New concept - pheno-match category. We only retain variants which are

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

    # for these categories, require a phenotype-gene match
    cats_require_pheno_match = config_retrieve(['category_rules', 'phenotype_match'], [])

    panel_meta: dict[int, str] = {content.id: content.name for content in panelapp_data.metadata}

    gene_details: dict[str, set[int]] = {}

    for each_event in result_list:
        # shouldn't be possible, here as a precaution
        assert each_event.categories, f'No categories for {each_event.var_data.coordinates.string_format}'

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
            phenotype_intersection = participant_panels.samples[each_event.sample].panels.intersection(all_panels)

            # is this gene relevant for this participant?
            # this test includes matched, cohort-level, and core panel
            if not phenotype_intersection.union(cohort_intersection):
                continue

            matched_panels = {
                panel_meta[pid]
                for pid in phenotype_intersection
                if pid != config_retrieve(['workflow', 'default_panel'], 137)
            }

        forced_panels = set()
        if cohort_intersection:
            forced_panels = {panel_meta[pid] for pid in cohort_intersection}

        # this is a horrible operation
        # if the variant-gene doesn't have a cohort-forced or phenotypic match panel
        # AND there's only one category assigned
        # AND that category is in the list of categories which require a phenotype match
        # skip this variant
        if (
            (not (forced_panels or matched_panels))
            and (len(each_event.support_vars) == 0)
            and (all(cat in cats_require_pheno_match for cat in each_event.categories))
        ):
            continue

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
            if prev_event.reasons != {'Autosomal Dominant'} and AMBIGUOUS_FLAG in both_flags:
                both_flags.remove(AMBIGUOUS_FLAG)
            prev_event.flags = both_flags

    # organise the variants by chromosomal location... why?
    for sample in results_holder.results:
        results_holder.results[sample].variants.sort()

    return results_holder


def count_families(pedigree: Pedigree, samples: set[str]) -> dict:
    """
    add metadata to results
    parsed during generation of the report

    affected, male, female, and family sizes all at the same level
    maybe re-think this output structure for the report

    Args:
        pedigree (Pedigree): the Pedigree object for the family
        samples (list): all the samples across all VCFs

    Returns:
        A breakdown of all the family structures within this analysis
    """

    trios_in_a_quad = 2

    # the final dict of counts to return
    stat_counter: dict[str, int] = defaultdict(int)

    # now count family sizes, structures, sexes, and affected
    for family_members in pedigree.by_family.values():
        # reduce those to participants within VCFs
        filter_members = [member for member in family_members if member.id in samples]

        # find trios. there's magic method for this in peds, but we're
        # using a different representation
        trios = [
            member
            for member in filter_members
            if member.affected == '2' and member.mother is not None and member.father is not None
        ]

        # this could be extended, or do more stringent family tests
        if len(trios) == 1:
            stat_counter['trios'] += 1
        elif len(trios) == trios_in_a_quad:
            stat_counter['quads'] += 1
        else:
            stat_counter[str(len(family_members))] += 1

        for member in filter_members:
            stat_counter[MALE_FEMALE[member.sex]] += 1
            if member.affected == '2':
                stat_counter['affected'] += 1

    return dict(stat_counter)


def prepare_results_shell(
    results_meta: ResultMeta,
    small_samples: set[str],
    sv_samples: set[str],
    pedigree: Pedigree,
    dataset: str,
    panelapp: PanelApp,
    panel_data: PhenotypeMatchedPanels | None = None,
) -> ResultData:
    """
    Creates a ResultData object, with participant metadata filled out

    Args:
        results_meta (): metadata for the results
        small_samples (): samples in the Small VCF
        sv_samples (): samples in the SV VCFs
        pedigree (): the Pedigree object
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

    # all affected samples in Pedigree, small variant and SV VCFs may not completely overlap
    all_samples: set[str] = small_samples | sv_samples

    for sample in [sam for sam in pedigree.members if sam.affected == '2' and sam.id in all_samples]:
        family_members = {
            member.id: FamilyMembers(
                **{'sex': MALE_FEMALE[member.sex], 'affected': member.affected == '2', 'ext_id': member.ext_id},
            )
            for member in pedigree.by_family[sample.family]
        }
        sample_panel_data = panel_data.samples.get(sample.id, ParticipantHPOPanels())
        results_shell.results[sample.id] = ParticipantResults(
            **{
                'variants': [],
                'metadata': ParticipantMeta(
                    **{
                        'ext_id': sample.ext_id,
                        'family_id': sample.family,
                        'members': family_members,
                        'phenotypes': sample_panel_data.hpo_terms,
                        'panel_ids': sample_panel_data.panels,
                        'panel_names': [panel_meta[panel_id] for panel_id in sample_panel_data.panels],
                        'solved': bool(sample.id in solved_cases or sample.family in solved_cases),
                        'present_in_small': sample.id in small_samples,
                        'present_in_sv': sample.id in sv_samples,
                    },
                ),
            },
        )

    return results_shell


def main(
    labelled_vcf: str,
    out_json: str,
    panelapp: str,
    pedigree: str,
    input_path: str,
    labelled_sv: list[str] | None = None,
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

    if dataset is None:
        dataset = config_retrieve(['workflow', 'dataset'])
    assert isinstance(dataset, str)

    if labelled_sv is None:
        labelled_sv = []

    out_json_path = to_path(out_json)

    # parse the pedigree from the file
    ped = make_flexible_pedigree(pedigree)

    # parse panelapp data from dict
    panelapp_data: PanelApp = read_json_from_path(panelapp, return_model=PanelApp)  # type: ignore

    # set up the inheritance checks
    moi_lookup = set_up_moi_filters(panelapp_data=panelapp_data, pedigree=ped)

    pheno_panels: PhenotypeMatchedPanels | None = read_json_from_path(
        participant_panels,
        return_model=PhenotypeMatchedPanels,  # type: ignore
        default=None,
    )

    # create the new gene map
    new_gene_map = get_new_gene_map(panelapp_data, pheno_panels, dataset)

    result_list: list[ReportVariant] = []

    # collect all sample IDs from each VCF type
    small_vcf_samples: set[str] = set()
    sv_vcf_samples: set[str] = set()

    # open the small variant VCF using a cyvcf2 reader
    vcf_opened = VCFReader(labelled_vcf)
    small_vcf_samples.update(set(vcf_opened.samples))

    # optional SV behaviour
    sv_opened = [VCFReader(sv_vcf) for sv_vcf in labelled_sv]
    for sv_vcf in sv_opened:
        sv_vcf_samples.update(set(sv_vcf.samples))

    all_samples: set[str] = small_vcf_samples.union(sv_vcf_samples)

    # obtain a set of all contigs with variants
    for contig in canonical_contigs_from_vcf(vcf_opened):
        # assemble {gene: [var1, var2, ..]}
        contig_dict = gather_gene_dict_from_contig(
            contig=contig,
            variant_source=vcf_opened,
            sv_sources=sv_opened,
            new_gene_map=new_gene_map,
            singletons=bool('singleton' in pedigree),
        )

        result_list.extend(
            apply_moi_to_variants(
                variant_dict=contig_dict,
                moi_lookup=moi_lookup,
                panelapp_data=panelapp_data.genes,
                pedigree=ped,
            ),
        )

    # do we have seqr projects?
    seqr_project = get_cohort_seq_type_conf().get('seqr_project')

    # create the full final output file
    results_meta = ResultMeta(
        **{
            'input_file': input_path,
            'cohort': dataset or config_retrieve(['workflow', 'dataset'], 'unknown'),
            'family_breakdown': count_families(ped, samples=all_samples),
            'panels': panelapp_data.metadata,
            'container': config_retrieve(['workflow', 'driver_image']),
            'projects': [seqr_project] if seqr_project else [],
        },
    )

    # create a shell to store results in, adds participant metadata
    results_model = prepare_results_shell(
        results_meta=results_meta,
        small_samples=small_vcf_samples,
        sv_samples=sv_vcf_samples,
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
    filter_results(results_model, singletons=bool('singleton' in pedigree), dataset=dataset)

    # write the output to long term storage using Pydantic
    # validate the model against the schema, then write the result if successful
    with to_path(out_json_path).open('w') as out_file:
        out_file.write(ResultData.model_validate(results_model).model_dump_json(indent=4))


if __name__ == '__main__':
    get_logger(__file__).info('Starting MOI testing phase')
    get_logger().info(f'Operational Categories: {CATEGORY_DICT}')

    parser = ArgumentParser(description='Startup commands for the MOI testing phase of Talos')
    parser.add_argument('--labelled_vcf', help='Category-labelled VCF')
    parser.add_argument('--labelled_sv', help='Category-labelled SV VCF', default=[], nargs='+')
    parser.add_argument('--out_json', help='Prefix to write JSON results to')
    parser.add_argument('--panelapp', help='Path to JSON file of PanelApp data')
    parser.add_argument('--pedigree', help='Path to joint-call PED file')
    parser.add_argument('--input_path', help='source data', default='Not supplied')
    parser.add_argument('--participant_panels', help='panels per participant', default=None)
    parser.add_argument('--dataset', help='optional, dataset to use', default=None)
    args = parser.parse_args()

    main(
        labelled_vcf=args.labelled_vcf,
        out_json=args.out_json,
        panelapp=args.panelapp,
        pedigree=args.pedigree,
        input_path=args.input_path,
        labelled_sv=args.labelled_sv,
        participant_panels=args.participant_panels,
        dataset=args.dataset,
    )
