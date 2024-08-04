"""
Takes the PanelApp data (the full region of interest)
produce an output of the {gene ID: gene symbol} for all included genes
write this to a file

This has been extracted to a separate script to allow for parallelisation
and to ensure that we can separately repeat this one step in isolation
if we run into API/throttling issues

This may  be integrated into the HPO~Phenotype matching script if it runs consistently enough
"""

import asyncio
import json
from argparse import ArgumentParser

from aiohttp import ClientSession

from talos.config import config_retrieve
from talos.models import PanelApp
from talos.utils import chunks, read_json_from_path

ENSEMBL_REST_API = 'http://rest.ensembl.org'


async def match_ensgs_to_symbols(genes: list[str], session: ClientSession) -> dict[str, str]:
    data_payload = json.dumps({'ids': genes})
    r = await session.request(
        method='POST',
        url=f'{ENSEMBL_REST_API}/lookup/id',
        headers={'Accept': 'application/json', 'Content-Type': 'application/json'},
        data=data_payload,
    )
    r.raise_for_status()
    json_reponse = await r.json()
    # match symbol to the ENSG (or Unknown if the key is missing, or has a None value)
    return {value.get('display_name'): key for key, value in json_reponse.items() if value}


async def match_symbol_to_ensg(gene_symbol: str, session: ClientSession) -> tuple[str, str]:
    r = await session.request(
        method='GET',
        url=f'{ENSEMBL_REST_API}/lookup/id/{gene_symbol}',
        headers={'Content-Type': 'application/json'},
    )
    r.raise_for_status()
    json_reponse = await r.json()
    return gene_symbol, json_reponse['display_name']


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--panelapp', help='Path to the PanelApp results file', required=True)
    parser.add_argument('--out_path', help='where to write the output (.json)', required=True)
    args = parser.parse_args()

    main(args.panelapp, out_path=args.out_path)


def main(panelapp_path: str, out_path: str):
    """

    Args:
        panelapp_path (str): path to the PanelApp model JSON
        out_path ():
    """
    panelapp_object = read_json_from_path(read_path=panelapp_path, return_model=PanelApp)
    # confirm for mypy
    assert isinstance(panelapp_object, PanelApp)
    genes = list(panelapp_object.genes.keys())
    ensg_to_symbol_mapping = asyncio.run(get_symbols_for_ensgs(genes))

    with open(out_path, 'w') as f:
        json.dump(ensg_to_symbol_mapping, f, indent=4)


async def get_symbols_for_ensgs(genes: list[str]) -> dict[str, str]:
    chunksize = config_retrieve(['FindGeneSymbolMap', 'chunk_size'], 800)
    all_results: dict[str, str] = {}
    async with ClientSession() as session:
        for chunk in chunks(genes, chunksize):
            all_results |= await match_ensgs_to_symbols(chunk, session=session)
    return all_results


if __name__ == '__main__':
    cli_main()
