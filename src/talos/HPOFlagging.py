"""
Takes the result date from ValidateMOI
For each reportable result, identify if the variant gene is enriched relateive to the participant's phenotypes
This uses the HPO3/pyHPO gene enrichment model, see:

https://github.com/anergictcell/pyhpo?tab=readme-ov-file#get-genes-enriched-in-an-hposet
https://pyhpo.readthedocs.io/en/latest/stats.html#statistics

This pretty cleanly matches our use-case: "You have a set of HPOTerms and want to find the most likely causative gene"

Here we're doing an enrichment test using the participant's phenotypes, then accepting any genes with an enrichment
probability below 0.05

1. Write 2 results - the result set annotated with Phenotype-match flags
2. The annotated result set, reduced to only phenotype matches

The removal of variants where a phenotype match is required but not found it now done here
A list of categories in config is used to determine which variants are required to have a phenotype match
Variants where all categories are on that list will be removed unless any of these are satisfied:
- variant gene is on a forced panel
- variant gene is on a phenotype match panel
- variant gene is phenotype matched to the family
"""

from argparse import ArgumentParser

from pyhpo import Gene, HPOSet, Ontology
from pyhpo.stats import EnrichmentModel

from talos.config import config_retrieve
from talos.models import ParticipantResults, ResultData
from talos.static_values import get_granular_date
from talos.utils import get_logger, phenotype_label_history, read_json_from_path


# some setup for hpo3
_ = Ontology()
ENRICHMENT_MODEL = EnrichmentModel('gene')


def get_gene_enrichment(hpo_terms: list[str]) -> dict[str, dict[str, float]]:
    """
    get the gene enrichment for a set of HPO terms
    Args:
        hpo_terms ():

    Returns:
        dictionary of genes and their enrichments
    """

    # pull the enrichment threshold from config
    enrichment_threshold = config_retrieve(['HPOFlagging', 'enrichment_threshold'], 0.05)

    # check for genes enriched in the participant's HPO terms
    participant_hpo3 = HPOSet.from_queries(hpo_terms)

    # find genes which are enriched relative to the participant's HPO terms
    filtered_enriched: dict[str, dict[str, float]] = {
        gene['item'].name: {'enrichment': gene['enrichment'], 'fold': gene['fold']}
        for gene in ENRICHMENT_MODEL.enrichment(method='hypergeom', hposet=participant_hpo3)
        if gene['enrichment'] <= enrichment_threshold
    }

    return filtered_enriched


def annotate_phenotype_matches(result_object: ResultData, ensg2symbol: str):
    """
    for each variant, find probability of the variant gene showing an enriched association with the
    set of patient phenotype terms

    Args:
        result_object (ResultData):
        ensg2symbol (str): path to the ENSG to symbol mapping
    """

    ensg_map: dict[str, str] = read_json_from_path(ensg2symbol)

    for participant in result_object.results.values():
        # no HPOs, no problems
        if not participant.metadata.phenotypes:
            continue

        participant_hpos_dict = {hpo.id: hpo.label for hpo in participant.metadata.phenotypes}
        participant_hpos = list(participant_hpos_dict.keys())

        # find genes which are enriched relative to the participant's HPO terms
        filtered_enriched = get_gene_enrichment(participant_hpos)

        for variant in participant.variants:
            var_symbol = ensg_map.get(variant.var_data.info['gene_id'])  # type: ignore[assignment]

            if var_symbol in filtered_enriched:
                variant.phenotype_labels.add(f'p.{filtered_enriched[var_symbol]["enrichment"]:.4f}')
                if variant.date_of_phenotype_match is None:
                    variant.date_of_phenotype_match = get_granular_date()

            # get exactly matching  HPO terms for the gene
            try:
                gene_hpo_set = Gene.get(var_symbol).hpo_set()
                gene_hpo_strings = {hpo_term.id for hpo_term in gene_hpo_set}
                for hpo_term in gene_hpo_strings & set(participant_hpos):
                    if variant.date_of_phenotype_match is None:
                        variant.date_of_phenotype_match = get_granular_date()
                    hpo_object = Ontology.get_hpo_object(hpo_term)
                    variant.phenotype_labels.add(f'{hpo_term}: {hpo_object.name} (exact)')
            except KeyError:
                get_logger(__file__).warning(f'Could not find gene {var_symbol} in the ontology')

    phenotype_label_history(result_object)


def remove_phenotype_required_variants(result_object: ResultData):
    """
    remove any variants where a phenotype match is required but not found

    Args:
        result_object ():
    """

    # for these categories, require a phenotype-gene match
    cats_require_pheno_match = config_retrieve(['ValidateMOI', 'phenotype_match'], [])

    for participant in result_object.results.values():
        kept_variants = []
        for variant in participant.variants:
            # boolean for whether this variant was phenotype matched
            matched_variant = bool(variant.phenotype_labels or variant.panels.matched or variant.panels.forced)

            # if the variant-gene doesn't have a cohort-forced or phenotypic match panel
            # and isn't a specific match between participant and gene symbol
            # AND all categories assigned required a phenotype match
            # AND the variant isn't a compound het
            # skip this variant
            if (
                (not matched_variant)
                and (len(variant.support_vars) == 0)
                and (all(cat in cats_require_pheno_match for cat in variant.categories))
            ):
                continue
            kept_variants.append(variant)

        participant.variants = kept_variants


def filter_and_write_out(annotated_results: ResultData, out_path: str):
    """
    filters the results in the Result model to only those with
    phenotypically matched variants

    Should update the metadata here, specifically the counting

    Args:
        annotated_results (ResultData):
        out_path (str):
    """

    new_rd = ResultData(version=annotated_results.version, metadata=annotated_results.metadata)

    for sample_id, sample_data in annotated_results.results.items():
        sample_participant_results = ParticipantResults(metadata=sample_data.metadata)

        for variant in sample_data.variants:
            if variant.phenotype_labels:
                sample_participant_results.variants.append(variant)

        if sample_participant_results.variants:
            new_rd.results[sample_id] = sample_participant_results

        # validate the object
        validated_results = ResultData.model_validate(new_rd)

        # validate and write using pydantic
        with open(out_path, 'w', encoding='utf-8') as handle:
            handle.write(validated_results.model_dump_json(indent=4))


def cli_main():
    parser = ArgumentParser(description='')
    parser.add_argument('--input', help='The Result data in JSON form')
    parser.add_argument('--gene_map', help='A map of gene symbol to ENSG for all genes in this analysis')
    parser.add_argument('--output', help='Annotated full output')
    parser.add_argument('--phenout', help='Annotated phenotype-only output')
    args = parser.parse_args()
    main(
        result_file=args.input,
        gene_map=args.gene_map,
        out_path=args.output,
        phenout=args.phenout,
    )


def main(
    result_file: str,
    gene_map: str,
    out_path: str,
    phenout: str | None = None,
):
    """

    Args:
        result_file (str): path to the ValidateMOI output JSON
        gene_map (str): output of FindGeneSymbolMap, JSON
        out_path (str): path to write the annotated results
        phenout (str): optional, path to phenotype filtered outputs
    """

    # read the results JSON into an object
    results = read_json_from_path(result_file, return_model=ResultData)

    # label phenotype matches, and write the results to a JSON file
    annotate_phenotype_matches(results, gene_map)

    # remove any variants where a phenotype match is required but not found
    remove_phenotype_required_variants(results)

    # validate the object
    validated_results = ResultData.model_validate(results)

    # validate and write using pydantic
    with open(out_path, 'w', encoding='utf-8') as handle:
        handle.write(validated_results.model_dump_json(indent=4))

    # reduce the JSON down to just phenotype matched variants, and the samples they occur in
    if phenout:
        filter_and_write_out(validated_results, out_path=phenout)


if __name__ == '__main__':
    cli_main()
