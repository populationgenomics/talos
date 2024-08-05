"""
Takes the result date from ValidateMOI
For each reportable result, identify if the variant is a strong phenotype match for the family

1. Write 2 results - the result set annotated with Phenotype-match flags
2. The annotated result set, reduced to only phenotype matches

Include something here about where to find the inputs ...

gene_to_phtnotype: Download: https://hpo.jax.org/app/data/annotations#:~:text=download-,GENES,-TO%20PHENOTYPE"
phenio.db: Download: https://data.monarchinitiative.org/monarch-kg/latest/phenio.db.gz
"""

from argparse import ArgumentParser
from collections import defaultdict

from semsimian import Semsimian

from talos.config import config_retrieve
from talos.models import ParticipantResults, ResultData
from talos.static_values import get_granular_date
from talos.utils import phenotype_label_history, read_json_from_path

_SEMSIM_CLIENT: Semsimian | None = None


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


def parse_genes_to_hpo(g2p_file: str, ensg2symbol: str, genes: set[str]) -> dict[str, set[str]]:
    """
    Parse genes to phenotype file from Jax.
    Returns a dict of gene_symbol -> set of HPO ids

    Args:
        g2p_file ():
        ensg2symbol (): the gene symbols to Ensembl gene IDs for genes relevant to this analysis
        genes (set[str]): all genes in this report

    Returns:
        dict[str, set[str]]: ENSG -> set of HPO ids
    """

    # lookup of {symbol: ENSG}
    ensg_map: dict[str, str] = read_json_from_path(ensg2symbol)

    gene_to_phenotype: dict[str, set[str]] = defaultdict(set)
    with open(g2p_file, encoding='utf-8') as f:
        for line in f:
            ncbi_gene_id, gene_symbol, hpo_id, hpo_name, frequency, disease_id = line.split('\t')
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
            ensgs.add(variant.var_data.info['gene_id'])

            # for structural variants, add all the LOF'd genes
            # TODO (mwelland): if/when we create other SV categories, we may need to catch those here too
            if lof := variant.var_data.info.get('predicted_lof'):
                ensgs.update(set(lof.split(',')))

    return ensgs


def annotate_phenotype_matches(result_object: ResultData, gen_phen: dict[str, set[str]], out_path: str):
    """
    for each variant, find any phenotype matches between the participant and gene HPO sets

    checks config for HPOFlagging.strict (bool)
    if strict, this will test for an exact overlap
    if strict is False, this will be done semantically via Semsimian

    Args:
        results (ResultData):
        gen_phen (dict): mapping of ENSGs to relevant HPO terms
        out_path (str): path to write results to
    """

    semantic_match = config_retrieve(['HPOFlagging', 'semantic_match'], False)

    min_similarity: float = config_retrieve(['HPOFlagging', 'min_similarity'])

    for participant in result_object.results.values():
        participant_hpos_dict = {hpo.id: hpo.label for hpo in participant.metadata.phenotypes}
        participant_hpos = set(participant_hpos_dict.keys())

        for variant in participant.variants:
            var_gene = variant.var_data.info['gene_id']
            gene_hpos = gen_phen.get(var_gene, set())

            # under strict matching we require exact overlapping terms
            # we always run a strict match
            for hpo_id in participant_hpos & gene_hpos:
                variant.phenotype_labels.add(f'{hpo_id}: {participant_hpos_dict[hpo_id]}')

            # optionally also use semantic matching for phenotypic similarity
            if participant_hpos and gene_hpos and semantic_match:
                termset_similarity = get_sem_client().termset_pairwise_similarity(participant_hpos, gene_hpos)
                # Convert object terms (gene_phenotypes) to lookup dict
                object_termset = {
                    term_dict['id']: term_dict['label']
                    for term in termset_similarity['object_termset']
                    for term_dict in term.values()
                }

                # Find all phenotype matches that meet the min_score threshold
                pheno_matches = {
                    f"{match['object_id']}: {object_termset[match['object_id']]}"
                    for match in termset_similarity['subject_best_matches']['similarity'].values()
                    if float(match['ancestor_information_content']) > min_similarity
                }
                if variant.date_of_phenotype_match is None:
                    variant.date_of_phenotype_match = get_granular_date()

                variant.phenotype_labels = pheno_matches

    phenotype_label_history(result_object)
    # validate the object
    validated_results = ResultData.model_validate(result_object)

    # validate and write using pydantic
    with open(out_path, 'w', encoding='utf-8') as handle:
        handle.write(validated_results.model_dump_json(indent=4))


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
    parser.add_argument('--results', help='The Result data in JSON form')
    parser.add_argument('--gene_map', help='A map of gene symbol to ENSG for all genes in this analysis')
    parser.add_argument('--gen2phen', help='path to the genotype-phenotype file')
    parser.add_argument('--phenio', help='A phenio DB file')
    parser.add_argument('--out', help='Annotated full output')
    parser.add_argument('--phenout', help='Annotated phenotype-only output')
    args = parser.parse_args()
    main(
        results=args.results,
        gene_map=args.gene_map,
        gen2phen=args.gen2phen,
        phenio=args.phenio,
        out_path=args.out,
        phenout=args.phenout,
    )


def main(
    results: str,
    gene_map: str,
    gen2phen: str,
    phenio: str,
    out_path: str,
    phenout: str | None = None,
):
    """

    Args:
        results (str): path to the ValidateMOI output JSON
        gene_map (str): output of FindGeneSymbolMap, JSON
        gen2phen (str): path to a test file of known Phenotypes per gene
        phenio (str): path to a PhenoIO DB file
        out_path (str): path to write the annotated results
        phenout (str): optional, path to phenotype filtered outputs
    """

    # read the results JSON into an object
    results = read_json_from_path(results, return_model=ResultData)

    # find all the genes present in this report
    relevant_ensgs = find_genes_in_these_results(results)

    # for each gene, identify the HPO terms that are associated with it
    gen_phen_dict = parse_genes_to_hpo(gen2phen, gene_map, genes=relevant_ensgs)

    # unless strict, set up the client with a phenio file
    # this is a bit hacky, but means we don't need to pass the phenio path around
    if not config_retrieve(['HPOFlagging', 'strict'], False):
        _semsim = get_sem_client(phenio_db=phenio)

    # label phenotype matches, and write the results to a JSON file
    annotate_phenotype_matches(results, gen_phen_dict, out_path=out_path)

    # reduce the JSON down to just phenotype matched variants, and the samples they occur in
    if phenout:
        filter_and_write_out(results, out_path=phenout)


if __name__ == '__main__':
    cli_main()
