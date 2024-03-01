"""
a script to post-process the AIP output data

produces a minimised representation of the output data,
containing only the data required for the SEQR app.

- Metadata, containing the description for each flag
- Results, containing per-individual results
    - Individual ID
        - Variant ID
            - Categories (list)
            - Support Variants (list)

Also produce a second version of the same, limited to phenotype-matches
"""

import json
import logging
from argparse import ArgumentParser

from reanalysis.models import MiniForSeqr, MiniVariant, ResultData


def coord_to_string(coord: dict) -> str:
    """
    converts a coordinate dict to a string

    Args:
        coord (dict): a coordinate dict

    Returns:
        str: a string representation of the coordinate dict
    """
    return f"{coord['chrom']}-{coord['pos']}-{coord['ref']}-{coord['alt']}"


def main(
    input_file: str, output: str, ext_map: str | None = None, pheno_match: bool = False
):
    """
    reads in the input file, shrinks it, and writes the output file


    Args:
        input_file (str):
        output (str):
        ext_map (str): optional mapping of internal to external IDs for seqr
        pheno_match (bool): whether to limit to phenotype-matching variants
    """

    with open(input_file, encoding='utf-8') as f:
        data = ResultData.model_validate(json.load(f))

    lil_data = MiniForSeqr(
        **{
            'metadata': {'categories': data.metadata.categories},
        }
    )
    ext_map_dict = None
    if ext_map:
        with open(ext_map, encoding='utf-8') as f:
            ext_map_dict = json.load(f)

    for individual, details in data.results.items():
        # optionally update to point to Seqr identities
        if ext_map_dict:
            individual = ext_map_dict.get(individual, individual)

        lil_data.results[individual] = {}
        for variant in details.variants:
            var_data = variant.var_data
            if pheno_match and not variant.panels.matched:
                continue
            lil_data.results[individual][var_data.info['seqr_link']] = MiniVariant(
                **{
                    'categories': variant.categories,
                    'support_vars': variant.support_vars,
                }
            )
    additional_string = 'phenotype-matched' if pheno_match else ''
    if not any(lil_data.results.values()):
        logging.info(f'No {additional_string} results found')
        return
    with open(output, 'w', encoding='utf-8') as f:
        f.write(MiniForSeqr.model_validate(lil_data).model_dump_json(indent=4))

    logging.info(f'Wrote {additional_string} output to {output}')


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('input_file', help='the input file to process')
    parser.add_argument('output_file', help='the output file to write to')
    parser.add_argument(
        '--external_map',
        help='mapping of internal to external IDs for seqr',
        default=None,
        type=str,
    )
    args = parser.parse_args()

    main(input_file=args.input_file, output=args.output_file, ext_map=args.external_map)
    main(
        input_file=args.input_file,
        output=args.output_file,
        ext_map=args.external_map,
        pheno_match=True,
    )
