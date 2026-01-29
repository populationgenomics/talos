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
from loguru import logger
from mendelbrot.pedigree_parser import PedigreeParser

from talos.config import config_retrieve
from talos.models import (
    FamilyMembers,
    MemberSex,
    PanelApp,
    PanelDetail,
    ParticipantHPOPanels,
    ParticipantMeta,
    ParticipantResults,
    ReportPanel,
    ReportVariant,
    ResultData,
    ResultMeta,
    translate_category,
)
from talos.moi_tests import MOIRunner
from talos.utils import (
    GeneDict,
    annotate_variant_dates_using_prior_results,
    canonical_contigs_from_vcf,
    find_comp_hets,
    gather_gene_dict_from_contig,
    generate_summary_stats,
    polish_exomiser_results,
    read_json_from_path,
)
from talos.version import __version__

AMBIGUOUS_FLAG = 'Ambiguous Cat.1 MOI'
MALE_FEMALE = {
    0: MemberSex.UNKNOWN.value,
    1: MemberSex.MALE.value,
    2: MemberSex.FEMALE.value,
}


def set_up_moi_filters(panelapp_data: PanelApp, pedigree: PedigreeParser) -> dict[str, MOIRunner]:
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
        pedigree (PedigreeParser):

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
    pedigree: PedigreeParser,
) -> list[ReportVariant]:
    """
    take all variants on a given contig & MOI filters
    find all variants/compound hets which fit the PanelApp MOI

    Args:
        variant_dict (dict): all possible variants, lists indexed by gene
        moi_lookup (dict): the MOI model runner per MOI string
        panelapp_data (dict): all genes and relevant details
        pedigree (PedigreeParser): the pedigree for this cohort
    """

    results = []

    for gene, variants in variant_dict.items():
        comp_het_dict = find_comp_hets(var_list=variants, pedigree=pedigree)

        # extract the panel data specific to this gene
        # extract once per gene, not once per variant
        panel_gene_data = panelapp_data.get(gene)

        # variant appears to be in a red gene
        if panel_gene_data is None:
            logger.error(f'How did this gene creep in? {gene}')
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
            if not isinstance(runner, MOIRunner):
                raise TypeError(f'MOIRunner was not a MOIRunner object: {runner}')

            variant_results = runner.run(
                principal_var=variant,
                comp_het=comp_het_dict,
                partial_pen=bool(variant.info.get('categorybooleanclinvarplp', False)),
            )

            # Flag! If this is a ClinVar P/LP variant, and we interpret under a lenient MOI, add flag for analysts
            if panel_gene_data.moi == 'Mono_And_Biallelic' and variant.info.get('categorybooleanclinvarplp', False):
                # consider each variant in turn
                for each_result in variant_results:
                    # never tag if this variant/sample is de novo
                    if translate_category('4') in each_result.categories:
                        continue

                    if each_result.reasons == 'Autosomal Dominant':
                        each_result.flags.add(AMBIGUOUS_FLAG)

            results.extend(variant_results)

    return results


def filter_results_to_panels(
    results_holder: ResultData,
    result_list: list[ReportVariant],
    panelapp: PanelApp,
) -> None:
    """
    It's possible 1 variant can be classified multiple ways
    e.g. different MOIs (dominant and comp het)
    e.g. the same variant with annotation from two genes

    Here we keep these all in parallel - some but not all annotations would be the same (locus, gnomAD, but not TX csq)

    stores names of panels specifically matched between the variant gene and participant:
    - as matched (phenotype matched)
    - as forced (cohort-wide applied panel)/custom panel

    Args:
        results_holder (): container for all results data
        result_list (): list of all ReportVariant events
        panelapp (PanelApp):

    Returns:
        None, object updated in-place
    """

    default_panel = config_retrieve(['GeneratePanelData', 'default_panel'], 137)

    # get the forced panel IDs from config, and add 0, the 'custom panel'
    forced_panel_ids = set(config_retrieve(['GeneratePanelData', 'forced_panels'], [])) | {0}

    # iterate over each separate reportable event
    for each_event in result_list:
        # find all panels featuring this gene
        gene_panels = panelapp.genes[each_event.gene].panels

        # get all forced panels this gene intersects with
        forced_panels_for_this_gene: set[int] = forced_panel_ids.intersection(gene_panels)

        # all panels assigned to this participant
        if each_event.sample in panelapp.participants:
            participant_panel_ids = set(panelapp.participants[each_event.sample].panels)
            # get union of naturally and forcibly matched panels for this gene, ignoring the default panel
            natural_matches_for_this_gene = participant_panel_ids.intersection(gene_panels) - {default_panel}
        else:
            logger.warning(f'Participant {each_event.sample} not found in panelapp participants')
            natural_matches_for_this_gene = set()

        # if this gene is not on a forced or naturally matched panel for this participant, skip
        if not (forced_panels_for_this_gene or natural_matches_for_this_gene or (default_panel in gene_panels)):
            continue

        # don't remove variants here, we do that in the pheno-matching stage
        each_event.panels = ReportPanel(
            matched={pid: panelapp.metadata[pid].name for pid in natural_matches_for_this_gene},
            forced={pid: panelapp.metadata[pid].name for pid in forced_panels_for_this_gene},
        )

        # add this event to the list for this participant
        results_holder.results[each_event.sample].variants.append(each_event)


def count_families(pedigree: PedigreeParser) -> dict:
    """
    add metadata to results
    parsed during generation of the report

    affected, male, female, and family sizes all at the same level
    maybe re-think this output structure for the report

    Args:
        pedigree (PedigreeParser): the Pedigree object for the family

    Returns:
        A breakdown of all the family structures within this analysis
    """

    trios_in_a_quad = 2

    # the final dict of counts to return
    stat_counter: dict[str | MemberSex, int] = defaultdict(int)

    # now count family sizes, structures, sexes, and affected
    for family_members in pedigree.by_family.values():
        # find trios
        trios = [
            member
            for member in family_members
            if member.is_affected and member.mother_id is not None and member.father_id is not None
        ]

        # count family size, based only on samples in the variant data
        if len(trios) == 1:
            stat_counter['trios'] += 1
        elif len(trios) == trios_in_a_quad:
            stat_counter['quads'] += 1
        else:
            stat_counter[str(len(family_members))] += 1

        for member in family_members:
            # any value not specifically addressed here defaults to unknown
            stat_counter[MALE_FEMALE.get(member.sex, MemberSex.UNKNOWN.value)] += 1
            if member.is_affected:
                stat_counter['affected'] += 1

    return dict(stat_counter)


def prepare_results_shell(
    results_meta: ResultMeta,
    small_samples: set[str],
    sv_samples: set[str],
    pedigree: PedigreeParser,
    panelapp: PanelApp,
) -> ResultData:
    """
    Creates a ResultData object, with participant metadata filled out

    Args:
        results_meta (): metadata for the results
        small_samples (): samples in the Small VCF
        sv_samples (): samples in the SV VCFs
        pedigree (): the Pedigree object, already reduced to samples in the callset
        panelapp (): dictionary of gene data

    Returns:
        ResultData with sample metadata filled in
    """

    # create an empty dict for all the samples
    results_shell = ResultData(metadata=results_meta)

    # find the solved cases in this project
    solved_cases = config_retrieve(['ValidateMOI', 'solved_cases'], [])

    for participant in pedigree.participants.values():
        if not participant.is_affected:
            continue

        family_members = {
            member.sample_id: FamilyMembers(
                sex=MALE_FEMALE.get(member.sex, MemberSex.UNKNOWN.value),
                affected=member.is_affected,
            )
            for member in pedigree.by_family[participant.family_id]
        }
        sample_panel_data = panelapp.participants.get(participant.sample_id, ParticipantHPOPanels())
        results_shell.results[participant.sample_id] = ParticipantResults(
            variants=[],
            metadata=ParticipantMeta(
                family_id=sample_panel_data.family_id or participant.family_id,
                members=family_members,
                phenotypes=sample_panel_data.hpo_terms,
                panel_details={
                    panel_id: panelapp.metadata[panel_id]
                    for panel_id in sample_panel_data.panels
                    if panel_id in panelapp.metadata
                },
                solved=bool(participant.sample_id in solved_cases or participant.family_id in solved_cases),
                present_in_small=participant.sample_id in small_samples,
                present_in_sv=participant.sample_id in sv_samples,
            ),
        )

    return results_shell


def cli_main():
    parser = ArgumentParser(description='Startup commands for the MOI testing phase of Talos')
    parser.add_argument('--labelled_vcf', help='Category-labelled VCF')
    parser.add_argument('--labelled_sv', help='Category-labelled SV VCF', default=None)
    parser.add_argument('--labelled_mito', help='Category-labelled Mito VCF', default=None)
    parser.add_argument('--output', help='Prefix to write JSON results to', required=True)
    parser.add_argument('--panelapp', help='QueryPanelApp JSON', required=True)
    parser.add_argument('--pedigree', help='Path to PED file', required=True)
    parser.add_argument('--previous', help='Path to previous results', default=None)
    args = parser.parse_args()

    main(
        labelled_vcf=args.labelled_vcf,
        output=args.output,
        panelapp_path=args.panelapp,
        pedigree=args.pedigree,
        labelled_sv=args.labelled_sv,
        labelled_mito=args.labelled_mito,
        previous=args.previous,
    )


def main(
    labelled_vcf: str,
    output: str,
    panelapp_path: str,
    pedigree: str,
    labelled_sv: str | None = None,
    labelled_mito: str | None = None,
    previous: str | None = None,
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
        labelled_mito (str | None): optional Mitochondrial VCF
        output (str): location to write output file
        panelapp_path (str): location of PanelApp data JSON
        pedigree (str): location of PED file
        previous (str | None): location of previous results JSON, or None if first time/history not required
    """
    logger.info(
        r"""Welcome To
    ███████████   █████████   █████          ███████     █████████
    █   ███   █  ███     ███   ███         ███     ███  ███     ███
        ███      ███     ███   ███        ███       ███ ███
        ███      ███████████   ███        ███       ███  █████████
        ███      ███     ███   ███        ███       ███         ███
        ███      ███     ███   ███      █  ███     ███  ███     ███
       █████    █████   █████ ███████████    ███████     █████████
        """,
    )

    panelapp: PanelApp = read_json_from_path(
        panelapp_path,
        return_model=PanelApp,
    )

    logger.info(f'Attempting to read history from {previous}')
    previous_results: ResultData | None = read_json_from_path(
        previous,
        return_model=ResultData,
        default=None,
    )

    result_list: list[ReportVariant] = []

    # collect all sample IDs from each VCF type
    small_vcf_samples: set[str] = set()
    sv_vcf_samples: set[str] = set()
    _mito_vcf_samples: set[str] = set()

    # open the small variant VCF using a cyvcf2 reader
    vcf_opened = VCFReader(labelled_vcf)
    small_vcf_samples.update(set(vcf_opened.samples))

    # optional SV behaviour
    sv_opened = None
    if labelled_sv:
        sv_opened = VCFReader(labelled_sv)
        sv_vcf_samples = set(sv_opened.samples)

    all_samples = small_vcf_samples.union(sv_vcf_samples)

    mito_opened = None
    if labelled_mito:
        mito_opened = VCFReader(labelled_mito)
        mito_vcf_samples = set(mito_opened.samples)
        all_samples |= mito_vcf_samples

    # parse the pedigree from the file
    ped = PedigreeParser(pedigree)

    # slim down the pedigree to only samples we have in the pedigrees
    ped.set_participants(ped.strip_pedigree_to_samples(all_samples))

    # reduce cohort to singletons, if the config says so
    if config_retrieve('singletons', False):
        logger.info('Reducing pedigree to affected singletons only')
        ped.set_participants(ped.as_singletons())
        ped.set_participants(ped.get_affected_members())

    # set up the inheritance checks
    moi_lookup = set_up_moi_filters(panelapp_data=panelapp, pedigree=ped)

    # obtain a set of all contigs with variants
    for contig in canonical_contigs_from_vcf(vcf_opened):
        # assemble {gene: [var1, var2, ..]}
        contig_dict = gather_gene_dict_from_contig(
            contig=contig,
            variant_source=vcf_opened,
            sv_source=sv_opened,
        )

        result_list.extend(
            apply_moi_to_variants(
                variant_dict=contig_dict,
                moi_lookup=moi_lookup,
                panelapp_data=panelapp.genes,
                pedigree=ped,
            ),
        )

    if mito_opened:
        contig_dict = gather_gene_dict_from_contig('chrM', variant_source=mito_opened)
        result_list.extend(
            apply_moi_to_variants(
                variant_dict=contig_dict,
                moi_lookup=moi_lookup,
                panelapp_data=panelapp.genes,
                pedigree=ped,
            ),
        )

    # create the full final output file
    results_meta = ResultMeta(
        family_breakdown=count_families(ped),
        panels=panelapp.metadata,
        version=__version__,
    )

    # create a shell to store results in, adds participant metadata
    results_model = prepare_results_shell(
        results_meta=results_meta,
        small_samples=small_vcf_samples,
        sv_samples=sv_vcf_samples,
        pedigree=ped,
        panelapp=panelapp,
    )

    # remove duplicate and invalid variants
    filter_results_to_panels(results_model, result_list, panelapp)

    # need some extra filtering here to tidy up exomiser categorisation
    polish_exomiser_results(results_model)

    # annotate previously seen results by building on previous analysis results
    annotate_variant_dates_using_prior_results(results_model, previous_results)

    generate_summary_stats(results_model)

    # write the output to long term storage using Pydantic
    # validate the model against the schema, then write the result if successful
    with open(output, 'w') as out_file:
        out_file.write(ResultData.model_validate(results_model).model_dump_json(indent=4))


if __name__ == '__main__':
    cli_main()
