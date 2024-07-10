"""
take the PanelApp data (the full region of interest)
take the participant & phenotype data
produce an output of the gene IDs with a phenotype match, and the HPO terms they are connected to

later we can use this in two ways
- we can feed phenotype-matched genes into a lower barrier to entry category
- we can pass over the report, and if a variant falls in a matched gene, and the patient has the phenotype, prioritise
"""

from argparse import ArgumentParser


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--panelapp', help='Path to the PanelApp results file', required=True)
    parser.add_argument('--phenio', help='path to the phenio db file', required=True)
    parser.add_argument('--output', help='where to write the output (.json)', required=True)
    args = parser.parse_args()
    print(args.phenio)


#     main(args.panelapp, phenio_path=args.phenio, output_path=args.output)
#
#
# def main(panelapp_path: str, phenio_path: str, output_path: str):
#     """
#
#     Args:
#         panelapp_path (str): path to the PanelApp model JSON
#         phenio_path ():
#         output_path ():
#
#     Returns:
#
#     """


if __name__ == '__main__':
    cli_main()
