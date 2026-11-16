"""
Takes the result date from ValidateMOI.
For each reportable result, identify if the variant is a strong phenotype match for the family

A list of categories in config is used to determine which variants are required to have a phenotype match
Variants where all categories are on that list will be removed unless any of these are satisfied:
- variant gene is on a forced panel
- variant gene is on a phenotype match panel
- variant gene is phenotype matched to the family
"""

from argparse import ArgumentParser
from collections import defaultdict

from semsimian import Semsimian

from talos.config import config_retrieve
from talos.models import PanelApp, ResultData, translate_category
from talos.static_values import get_granular_date
from talos.utils import get_symbol_to_ensg_mapping, read_json_from_path

_SEMSIM_CLIENT: Semsimian | None = None

DEFAULT_SEMANTIC_MATCH = False
# as documented in docs/Configuration.md
DEFAULT_MIN_SIMILARITY = 14.0


def get_sem_client(phenio_db: str | None = None) -> Semsimian:
    """
    create or retrieve a Semsimian client

    Args:
        phenio_db (str | None): needs to be present for the first call

    Returns:
        Semsimian instance
    """
    global _SEMSIM_CLIENT
    if _SEMSIM_CLIENT is None:
        _SEMSIM_CLIENT = Semsimian(spo=None, predicates=['rdfs:subClassOf'], resource_path=phenio_db)
    return _SEMSIM_CLIENT


def parse_genes_to_hpo(g2p_file: str, ensg_map: dict[str, str], genes: set[str]) -> dict[str, set[str]]:
    """
    Parse genes to phenotype file from Jax.
    Returns a dict of gene_symbol -> set of HPO ids

    Args:
        g2p_file ():
        ensg_map (): gene symbols to Ensembl gene IDs
        genes (set[str]): all genes in this report

    Returns:
        dict[str, set[str]]: ENSG -> set of HPO ids
    """

    gene_to_phenotype: dict[str, set[str]] = defaultdict(set)
    with open(g2p_file, encoding='utf-8') as f:
        for line in f:
            _, gene_symbol, hpo_id, _, _, _ = line.split('\t')
            if gene_symbol in ensg_map and ensg_map[gene_symbol] in genes:
                gene_to_phenotype.setdefault(ensg_map[gene_symbol], set()).add(hpo_id)
    return gene_to_phenotype


def find_genes_in_these_results(result_object: ResultData) -> set[str]:
    """
    identify all the genes we care about finding phenotype matches for

    Args:
        result_object (ResultData):

    Returns:
        set[str], all ensembl gene IDs in
    """
    ensgs: set[str] = set()
    for participant in result_object.results.values():
        for variant in participant.variants:
            # get the gene ID
            ensgs.add(variant.var_data.info['gene_id'])  # type: ignore[arg-type]

            # for structural variants, add all the LOF'd genes
            if lof := variant.var_data.info.get('lof'):
                ensgs.update(set(str(lof).split(',')))

    return ensgs


def semantic_match_enabled() -> bool:
    """HPOFlagging.semantic_match - one lookup, so client setup and matching can't disagree on the default"""
    return config_retrieve(['HPOFlagging', 'semantic_match'], DEFAULT_SEMANTIC_MATCH)


def semantic_matches(participant_hpos: set[str], gene_hpos: set[str], min_similarity: float) -> set[str]:
    """
    gene HPO terms semantically similar to the participant's, via Semsimian

    Args:
        participant_hpos (set[str]): the participant's HPO IDs
        gene_hpos (set[str]): HPO IDs associated with the variant's gene
        min_similarity (float): ancestor information content a match must exceed

    Returns:
        set[str], labels as 'HP:ID: label'
    """
    termset_similarity = get_sem_client().termset_pairwise_similarity(participant_hpos, gene_hpos)
    # Convert object terms (gene_phenotypes) to lookup dict
    object_termset = {
        term_dict['id']: term_dict['label']
        for term in termset_similarity['object_termset']
        for term_dict in term.values()
    }

    # Find all phenotype matches that meet the min_score threshold
    return {
        f'{match["object_id"]}: {object_termset[match["object_id"]]}'
        for match in termset_similarity['subject_best_matches']['similarity'].values()
        if float(match['ancestor_information_content']) > min_similarity
    }


def annotate_phenotype_matches(result_object: ResultData, gen_phen: dict[str, set[str]]):
    """
    for each variant, find any phenotype matches between the participant and gene HPO sets

    An exact overlap of HPO terms always counts as a match. If HPOFlagging.semantic_match is set, semantically
    similar terms (via Semsimian) are added on top. The first time any match is found, the variant gets today's
    date as its date_of_phenotype_match - history is applied in ValidateMOI, before this stage, so a date already
    present is from an earlier cycle and is kept.

    Args:
        result_object (ResultData):
        gen_phen (dict): mapping of ENSGs to relevant HPO terms
    """
    semantic_match = semantic_match_enabled()

    # only consulted for semantic matching, so only required for it
    min_similarity: float = DEFAULT_MIN_SIMILARITY
    if semantic_match:
        min_similarity = config_retrieve(['HPOFlagging', 'min_similarity'], DEFAULT_MIN_SIMILARITY)

    for participant in result_object.results.values():
        participant_hpos_dict = {hpo.id: hpo.label for hpo in participant.metadata.phenotypes}
        participant_hpos = set(participant_hpos_dict.keys())

        for variant in participant.variants:
            var_gene = variant.var_data.info['gene_id']
            gene_hpos = gen_phen.get(var_gene, set())  # type: ignore[arg-type]

            # an exact overlap always counts
            matches = {f'{hpo_id}: {participant_hpos_dict[hpo_id]}' for hpo_id in participant_hpos & gene_hpos}

            # semantic matches are added on top of the exact ones, never in place of them
            if semantic_match and participant_hpos and gene_hpos:
                matches |= semantic_matches(participant_hpos, gene_hpos, min_similarity)

            # don't assign a date if there are no matches
            if not matches:
                continue

            variant.phenotype_labels |= matches

            if variant.date_of_phenotype_match is None:
                variant.date_of_phenotype_match = get_granular_date()


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
                and (all(translate_category(cat) in cats_require_pheno_match for cat in variant.categories))
            ):
                continue
            kept_variants.append(variant)

        participant.variants = kept_variants


def cli_main():
    parser = ArgumentParser(description='')
    parser.add_argument('--input', help='The Result data in JSON form')
    parser.add_argument('--panelapp', help='JSON of PanelApp data')
    parser.add_argument('--gen2phen', help='path to the genotype-phenotype file')
    parser.add_argument('--phenio', help='A phenio DB file')
    parser.add_argument('--output', help='Annotated output')
    parser.add_argument('--previous', help='Path to previous results', default=None)
    args = parser.parse_args()
    main(
        result_file=args.input,
        panelapp_path=args.panelapp,
        gen2phen=args.gen2phen,
        phenio=args.phenio,
        out_path=args.output,
    )


def main(
    result_file: str,
    panelapp_path: str,
    gen2phen: str,
    phenio: str,
    out_path: str,
):
    """

    Args:
        result_file (str): path to the ValidateMOI output JSON
        panelapp_path (str): JSON of PanelApp data in this analysis
        gen2phen (str): path to a test file of known Phenotypes per gene
        phenio (str): path to a PhenoIO DB file
        out_path (str): path to write the annotated results
    """

    panelapp = read_json_from_path(panelapp_path, return_model=PanelApp)

    gene_map = get_symbol_to_ensg_mapping(panelapp)

    # read the results JSON into an object
    results = read_json_from_path(result_file, return_model=ResultData)

    # find all the genes present in this report
    relevant_ensgs = find_genes_in_these_results(results)

    # for each gene, identify the HPO terms that are associated with it
    gen_phen_dict = parse_genes_to_hpo(gen2phen, gene_map, genes=relevant_ensgs)

    # only semantic matching needs the client - set it up with the phenio file here
    # this is a bit hacky, but means we don't need to pass the phenio path around
    if semantic_match_enabled():
        get_sem_client(phenio_db=phenio)

    # label phenotype matches, and write the results to a JSON file
    annotate_phenotype_matches(results, gen_phen_dict)

    # remove any variants where a phenotype match is required but not found
    remove_phenotype_required_variants(results)

    # validate the object
    validated_results = ResultData.model_validate(results)

    # validate and write using pydantic
    with open(out_path, 'w', encoding='utf-8') as handle:
        handle.write(validated_results.model_dump_json(indent=4))


if __name__ == '__main__':
    cli_main()
