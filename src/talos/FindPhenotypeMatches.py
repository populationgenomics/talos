# """
# take the PanelApp data (the full region of interest)
# take the participant & phenotype data
# produce an output of the gene IDs with a phenotype match, and the HPO terms they are connected to
#
# later we can use this in two ways
# - we can feed phenotype-matched genes into a lower barrier to entry category
# - we can pass over the report, and if a variant falls in a matched gene, and the patient has the phenotype, prioritise
# """
#
# from argparse import ArgumentParser
# from collections import defaultdict
#
# from semsimian import Semsimian
#
# from talos.models import PanelApp
# from talos.utils import read_json_from_path
#
#
# def cli_main():
#     parser = ArgumentParser()
#     parser.add_argument('--panelapp', help='Path to the PanelApp results file', required=True)
#     parser.add_argument('--g2p', help='path to the genotype-phenotype file', required=True)
#     parser.add_argument('--phenio', help='path to the phenio db file', required=True)
#     parser.add_argument('--output', help='where to write the output (.json)', required=True)
#     args = parser.parse_args()
#
#     main(args.panelapp, gen_to_phen_path=args.g2p, phenio_path=args.phenio, output_path=args.output)
#
#
# def parse_genes_to_phenotype(genes_to_phenotype_file: str) -> dict[str, set[str]]:
#     """
#     Parse genes to phenotype file from Jax.
#     Returns a dict of gene_symbol -> set of HPO ids
#     """
#     gene_to_phenotype = defaultdict(set)
#     with open(genes_to_phenotype_file) as f:
#         for line in f:
#             ncbi_gene_id, gene_symbol, hpo_id, hpo_name, frequency, disease_id = line.split('\t')
#             gene_to_phenotype[gene_symbol].add(hpo_id)
#     return gene_to_phenotype
#
#
# def main(panelapp_path: str, gen_to_phen_path: str, phenio_path: str, output_path: str, min_similarity: float = 14.0):
#     """
#
#     Args:
#         panelapp_path (str): path to the PanelApp model JSON
#         gen_to_phen_path (str): path to teh
#         phenio_path ():
#         output_path ():
#         min_similarity (float): min similarity to associate a gene and HPO term
#
#     Returns:
#
#     """
#     _panelapp_object = read_json_from_path(read_path=panelapp_path, return_model=PanelApp)
#     print(phenio_path)
#     sem_client = Semsimian(spo=None, predicates=['rdfs:subClassOf'], resource_path=phenio_path)
#     print(gen_to_phen_path)
#     gen_phen_dict = parse_genes_to_phenotype(gen_to_phen_path)
#     print(output_path)
#     print(min_similarity)
#
#
# if __name__ == '__main__':
#     cli_main()
