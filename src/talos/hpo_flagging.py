"""
Takes the result date from ValidateMOI.
For each reportable result, identify if the variant is a strong phenotype match for the family

Only terms within HPO's Phenotypic abnormality subtree take part in matching - modes of inheritance, onset and other
clinical modifiers say nothing about whether a gene fits a participant.

A list of categories in config is used to determine which variants are required to have a phenotype match
Variants where all categories are on that list will be removed unless any of these are satisfied:
- variant gene is on a forced panel
- variant gene is on a phenotype match panel
- variant gene is phenotype matched to the family
"""

from argparse import ArgumentParser
from collections import defaultdict

import networkx as nx
from loguru import logger
from obonet import read_obo

from talos.config import config_retrieve
from talos.models import PanelApp, ResultData, translate_category
from talos.static_values import get_granular_date
from talos.utils import get_symbol_to_ensg_mapping, read_json_from_path

PHENOTYPIC_ABNORMALITY = 'HP:0000118'

DEFAULT_SEMANTIC_MATCH = True
# as documented in docs/Configuration.md - equivalent to the Semsimian ancestor IC > 14.0 previously used against Phenio
DEFAULT_MAX_ANCESTOR_DESCENDANTS = 25


class HpoOntology:
    """
    The is_a hierarchy of HPO's Phenotypic abnormality subtree, as used for phenotype matching

    Semantic similarity is Resnik's: two terms are as similar as their most specific shared ancestor, and an
    ancestor's specificity is how few terms (itself included) sit beneath it.
    """

    def __init__(self, obo_path: str):
        graph = read_obo(obo_path, ignore_obsolete=False)

        self.names: dict[str, str] = {}
        # alt_ids and obsolete-with-replacement IDs, pointing at the current term
        self._redirect: dict[str, str] = {}
        for term, data in graph.nodes(data=True):
            self.names[term] = data.get('name', term)
            if data.get('is_obsolete') == 'true':
                if replaced_by := data.get('replaced_by'):
                    self._redirect[term] = replaced_by[0]
                continue
            for alt_id in data.get('alt_id', []):
                self._redirect[alt_id] = term

        # edges run child -> parent
        is_a = nx.DiGraph((child, parent) for child, parent, key in graph.edges(keys=True) if key == 'is_a')

        # ontology descendants of the root, which networkx calls ancestors, as edges point upwards
        subtree = nx.ancestors(is_a, PHENOTYPIC_ABNORMALITY) | {PHENOTYPIC_ABNORMALITY}
        is_a = is_a.subgraph(subtree)

        # every term's ancestors, itself included - parents are resolved before their children
        self.ancestors: dict[str, frozenset[str]] = {}
        for term in reversed(list(nx.topological_sort(is_a))):
            self.ancestors[term] = frozenset(
                {term}.union(*(self.ancestors[parent] for parent in is_a.successors(term)))
            )

        self.descendant_count: dict[str, int] = defaultdict(int)
        for term_ancestors in self.ancestors.values():
            for ancestor in term_ancestors:
                self.descendant_count[ancestor] += 1

    def phenotypic_abnormalities(self, terms: set[str]) -> set[str]:
        """
        current IDs for each term within Phenotypic abnormality (the root itself excluded) - others are dropped

        Args:
            terms (set[str]): HPO IDs, which may be alt_ids or obsolete IDs with a replacement

        Returns:
            set[str]
        """
        current = {self._redirect.get(term, term) for term in terms}
        return {term for term in current if term in self.ancestors and term != PHENOTYPIC_ABNORMALITY}

    def label(self, term: str) -> str:
        return f'{term}: {self.names.get(term, term)}'

    def ancestor_index(self, terms: set[str]) -> dict[str, str]:
        """
        every ancestor of these terms, mapped to the lowest-ID term beneath it - built once per gene

        Args:
            terms (set[str]): current, Phenotypic abnormality HPO IDs, e.g. a gene's annotations

        Returns:
            dict[str, str]: ancestor -> one of the terms it covers
        """
        index: dict[str, str] = {}
        for term in sorted(terms, reverse=True):
            for ancestor in self.ancestors[term]:
                index[ancestor] = term
        return index

    def best_match(self, term: str, index: dict[str, str], max_ancestor_descendants: int) -> str | None:
        """
        the indexed term most similar to this one, if their most specific shared ancestor is specific enough

        Args:
            term (str): current, Phenotypic abnormality HPO ID, e.g. a participant's
            index (dict[str, str]): from ancestor_index
            max_ancestor_descendants (int): the shared ancestor may have at most this many descendants, itself included

        Returns:
            str | None: the matched term, or None if there's no match
        """
        shared = [ancestor for ancestor in self.ancestors[term] if ancestor in index]
        # the root is always shared, so this is never empty
        best = min(shared, key=lambda ancestor: (self.descendant_count[ancestor], ancestor))
        if self.descendant_count[best] > max_ancestor_descendants:
            return None
        return index[best]


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


def annotate_phenotype_matches(result_object: ResultData, gen_phen: dict[str, set[str]], ontology: HpoOntology):
    """
    for each variant, find any phenotype matches between the participant and gene HPO sets

    Only Phenotypic abnormality terms are compared. An exact overlap of terms always counts as a match. If
    HPOFlagging.semantic_match is set, each participant term without an exact match is also matched to its most
    similar gene term, when their most specific shared ancestor has no more than
    HPOFlagging.max_ancestor_descendants descendants. The first time any match is found, the variant gets today's
    date as its date_of_phenotype_match - history is applied in ValidateMOI, before this stage, so a date already
    present is from an earlier cycle and is kept.

    Args:
        result_object (ResultData):
        gen_phen (dict): mapping of ENSGs to relevant HPO terms
        ontology (HpoOntology): the HPO hierarchy
    """
    semantic_match = config_retrieve(['HPOFlagging', 'semantic_match'], DEFAULT_SEMANTIC_MATCH)
    max_ancestor_descendants = config_retrieve(
        ['HPOFlagging', 'max_ancestor_descendants'],
        DEFAULT_MAX_ANCESTOR_DESCENDANTS,
    )
    if config_retrieve(['HPOFlagging', 'min_similarity'], None) is not None:
        logger.warning('HPOFlagging.min_similarity is no longer used, see HPOFlagging.max_ancestor_descendants')

    gene_terms = {gene: ontology.phenotypic_abnormalities(terms) for gene, terms in gen_phen.items()}
    gene_indexes: dict[str, dict[str, str]] = {}

    for participant in result_object.results.values():
        participant_terms = ontology.phenotypic_abnormalities({hpo.id for hpo in participant.metadata.phenotypes})

        for variant in participant.variants:
            var_gene: str = variant.var_data.info['gene_id']  # type: ignore[assignment]
            this_gene_terms = gene_terms.get(var_gene, set())

            # an exact overlap always counts
            matches = participant_terms & this_gene_terms

            # semantic matches are only sought for participant terms without an exact match
            if semantic_match and this_gene_terms:
                if var_gene not in gene_indexes:
                    gene_indexes[var_gene] = ontology.ancestor_index(this_gene_terms)
                for term in participant_terms - matches:
                    if match := ontology.best_match(term, gene_indexes[var_gene], max_ancestor_descendants):
                        matches.add(match)

            # don't assign a date if there are no matches
            if not matches:
                continue

            variant.phenotype_labels |= {ontology.label(term) for term in matches}

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
    parser.add_argument('--hpo', help='HPO ontology in OBO format', required=True)
    parser.add_argument('--output', help='Annotated output')
    parser.add_argument('--previous', help='Path to previous results', default=None)
    args = parser.parse_args()
    main(
        result_file=args.input,
        panelapp_path=args.panelapp,
        gen2phen=args.gen2phen,
        hpo_file=args.hpo,
        out_path=args.output,
    )


def main(
    result_file: str,
    panelapp_path: str,
    gen2phen: str,
    hpo_file: str,
    out_path: str,
):
    """

    Args:
        result_file (str): path to the ValidateMOI output JSON
        panelapp_path (str): JSON of PanelApp data in this analysis
        gen2phen (str): path to a test file of known Phenotypes per gene
        hpo_file (str): path to the HPO ontology OBO file, the same release as gen2phen
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

    # label phenotype matches
    annotate_phenotype_matches(results, gen_phen_dict, HpoOntology(hpo_file))

    # remove any variants where a phenotype match is required but not found
    remove_phenotype_required_variants(results)

    # validate the object
    validated_results = ResultData.model_validate(results)

    # validate and write using pydantic
    with open(out_path, 'w', encoding='utf-8') as handle:
        handle.write(validated_results.model_dump_json(indent=4))


if __name__ == '__main__':
    cli_main()
