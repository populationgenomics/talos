#!/usr/bin/env python3

"""
Not currently in use - this doesn't scale appropriately for our needs at this time

Take the following inputs:
    - a lookup of the gene symbol to the Ensembl gene ID for each gene in the ROI
    - a lookup of the gene symbol to phenotypic HPO terms
    ...
"""

from argparse import ArgumentParser
from collections import defaultdict

from semsimian import Semsimian

from talos.models import PhenotypeMatchedPanels
from talos.utils import read_json_from_path


def parse_genes_to_phenotype(g2p_file: str, symbols_to_ensembl: dict[str, str]) -> dict[frozenset, set[str]]:
    """
    Parse genes to phenotype file from Jax.
    Returns a dict of gene_symbol -> set of HPO ids

    Args:
        g2p_file ():
        symbols_to_ensembl (): the gene symbols to Ensembl gene IDs for genes relevant to this analysis

    Returns:
        dict[str, set[str]]: gene symbol -> set of HPO ids
    """

    gene_symbols_in_analysis = set(symbols_to_ensembl.values())
    gene_to_phenotype = defaultdict(set)
    with open(g2p_file) as f:
        for line in f:
            ncbi_gene_id, gene_symbol, hpo_id, hpo_name, frequency, disease_id = line.split('\t')
            if gene_symbol in gene_symbols_in_analysis:
                gene_to_phenotype[gene_symbol].add(hpo_id)

    # now do a hot rearrangement of the dict - we catch all unique sets of HPOs
    # in experiment this has reduced the search space from ~4k symbols to ~3k unique HPO sets
    results: dict[frozenset, set[str]] = defaultdict(set)
    for gene, hpos in gene_to_phenotype.items():
        results[frozenset(hpos)].add(gene)

    return results


def collect_unique_hpo_groups(party_panels: PhenotypeMatchedPanels) -> dict[frozenset, set[str]]:
    """
    Collect unique HPO groups from the ParticipantHPOPanels object
    The intention here is to reduce the number of individual comparisions being made
    Args:
        party_panels (PhenotypeMatchedPanels): the ParticipantHPOPanels object

    Returns:
        dict[frozenset, set[str]]: a dict of frozenset of HPO ids -> set of gene symbols
    """

    uniq_hpo_groups = defaultdict(set)
    for participant_id, content in party_panels.samples.items():
        if not content.hpo_terms:
            continue
        hpo_frozen = frozenset(term.id for term in content.hpo_terms)
        uniq_hpo_groups[hpo_frozen].add(participant_id)

    return uniq_hpo_groups


def find_similar_genes(
    gen_phen_dict: dict[frozenset, set[str]],
    hpo_groups,
    sem_client,
    min_similarity: float,
) -> dict[frozenset, dict[str, set[str]]]:
    """
    Find the genes that are similar to the HPO terms in the ParticipantHPOPanels object
    Args:
        gen_phen_dict (dict[str, set[str]]): gene symbol -> set of HPO ids
        hpo_groups (dict[frozenset, set[str]]): a dict of frozenset of HPO ids -> set of gene symbols
        sem_client (Semsimian): a SemSimian client
        min_similarity (float): min similarity to associate a gene and HPO term

    Returns:
        dict[frozenset, dict[str, set[str]]]: a dict of frozenset of HPO ids -> set of gene symbols which are similar
    """
    results: dict[frozenset, dict[str, set[str]]] = defaultdict(dict)
    for hpo_group in hpo_groups:
        for gene_hpos, gene_symbols in gen_phen_dict.items():
            termset_similarity = sem_client.termset_pairwise_similarity(set(hpo_group), set(gene_hpos))

            # Convert object terms (gene_phenotypes) to lookup dict
            object_termset = {
                term_dict['id']: term_dict['label']
                for term in termset_similarity['object_termset']
                for term_dict in term.values()
            }

            # Find all phenotype matches that meet the min_score threshold
            matches = {
                f"{match['object_id']}: {object_termset[match['object_id']]}"
                for match in termset_similarity['subject_best_matches']['similarity'].values()
                if float(match['ancestor_information_content']) > min_similarity
            }
            # this is a pretty non-specific way of handling this
            # we return all reasons we matched this collection of HPO terms, and the corresponding genes
            if matches:
                results[hpo_group].setdefault('symbols', set()).update(gene_symbols)
                results[hpo_group].setdefault('context', set()).update(matches)

    return results


def reassociate_participants_with_matched_genes(
    matched_hpo_to_genes: dict[frozenset, dict[str, set[str]]],
    hpo_groups: dict[frozenset, set[str]],
    party_panels: PhenotypeMatchedPanels,
):
    """
    Take the matched genes, and associate them with the participants
    Args:
        matched_hpo_to_genes (dict[frozenset, set[str]]): a dict of frozenset of HPO ids -> set of gene symbols
        hpo_groups (dict[frozenset, set[str]]): a dict of frozenset of HPO ids -> set of participant IDs
        party_panels (PhenotypeMatchedPanels): the ParticipantHPOPanels object
    """

    for frozen_hpo, content in matched_hpo_to_genes.items():
        for participant in hpo_groups[frozen_hpo]:
            party_panels.samples[participant].matched_genes.update(content['symbols'])
            party_panels.samples[participant].matched_phenotypes.update(content['context'])


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--g2p', help='path to the genotype-phenotype file', required=True)
    parser.add_argument('--gs2id', help='path to JSON with gene symbol to Ensembl gene ID', required=True)
    parser.add_argument('--panels', help='path to the participant-HPO-panels file', required=True)
    parser.add_argument('--phenio', help='path to the phenio db file', required=True)
    parser.add_argument('--output', help='where to write the output (.json)', required=True)
    args = parser.parse_args()
    main(gen2phen_path=args.g2p, gs2id=args.gs2id, panels_path=args.panels, phenio_db=args.phenio, out_path=args.output)


def main(gen2phen_path: str, gs2id: str, panels_path: str, phenio_db: str, out_path: str, min_similarity: float = 14.0):
    """

    Args:
        gen2phen_path (str): path to the
        gs2id (str): path to the JSON with gene symbol to Ensembl gene ID
        panels_path (str): path to the participant-HPO-panels file
        phenio_db ():
        out_path ():
        min_similarity (float): min similarity to associate a gene and HPO term
    """
    # read the JSON file containing the gene symbol to Ensembl gene ID
    symbols_to_ensembl = read_json_from_path(gs2id)

    # read the ParticipantHPOPanels object
    party_panels = read_json_from_path(panels_path, return_model=PhenotypeMatchedPanels)
    assert isinstance(party_panels, PhenotypeMatchedPanels)

    # identify all the unique HPO groups, and the samples that are in each group
    hpo_groups = collect_unique_hpo_groups(party_panels)

    # for each gene, identify the HPO terms that are associated with it
    gen_phen_dict = parse_genes_to_phenotype(gen2phen_path, symbols_to_ensembl)

    # create a SemSimian client
    sem_client = Semsimian(spo=None, predicates=['rdfs:subClassOf'], resource_path=phenio_db)

    # run the tests! This may be a horribly inefficient way to do this...
    matched_hpo_to_genes = find_similar_genes(gen_phen_dict, hpo_groups, sem_client, min_similarity)

    # populate the ParticipantHPOPanels object with the matched genes per participant, and reason for matches
    reassociate_participants_with_matched_genes(matched_hpo_to_genes, hpo_groups, party_panels)

    # validate the object
    valid_pheno_dict = PhenotypeMatchedPanels.model_validate(party_panels)

    # validate and write using pydantic
    with open(out_path, 'w', encoding='utf-8') as handle:
        handle.write(valid_pheno_dict.model_dump_json(indent=4))


if __name__ == '__main__':
    cli_main()
