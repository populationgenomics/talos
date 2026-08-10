#!/usr/bin/env python3

"""
Takes a mito VCF annotated with BCFtools csq and echtvar (MitImpact, MitoTip, nAPOGEE), applies fresh
ClinvArbitration decisions, and writes a labelled VCF containing only ClinVar P/LP variants in green mito genes.

This is a streaming (cyvcf2) process - no Hail, no Spark. One variant is held in memory at a time.
The ClinVar decisions are read from the ClinvArbitration TSV, restricted to mitochondrial contigs.
"""

import sys
from argparse import ArgumentParser

from cyvcf2 import VCF, Writer
from loguru import logger

from talos.models import PanelApp
from talos.utils import read_json_from_path
from talos.vcf_streaming import (
    PATHOGENIC,
    consequences_to_csq_string,
    normalise_chrom,
    parse_bcsq_entries,
    parse_pedigree,
    read_clinvar_decisions_tsv,
    split_csq_header,
    subset_reader_to_pedigree,
    variant_is_pass,
)

# mitochondrial contig names, in normalised (chr-less, MT-as-M) form
MITO_CONTIGS = {'M'}

# INFO fields added to the labelled output
OUTPUT_INFO_HEADERS = [
    {'ID': 'clinvar_significance', 'Number': '1', 'Type': 'String', 'Description': 'ClinvArbitration significance'},
    {'ID': 'clinvar_stars', 'Number': '1', 'Type': 'Integer', 'Description': 'ClinvArbitration gold stars'},
    {'ID': 'clinvar_allele', 'Number': '1', 'Type': 'Integer', 'Description': 'ClinvArbitration allele ID'},
    {'ID': 'clinvar_talos', 'Number': '1', 'Type': 'Integer', 'Description': 'ClinvArbitration Pathogenic flag'},
    {'ID': 'categorybooleanclinvarplp', 'Number': '1', 'Type': 'Integer', 'Description': 'ClinVar P/LP, 1+ stars'},
    {'ID': 'gnomad_AC', 'Number': '1', 'Type': 'Integer', 'Description': 'gnomAD AC (unpopulated for mito)'},
    {'ID': 'gnomad_AF', 'Number': '1', 'Type': 'Float', 'Description': 'gnomAD AF (unpopulated for mito)'},
    {'ID': 'gnomad_AC_XY', 'Number': '1', 'Type': 'Integer', 'Description': 'gnomAD AC XY (unpopulated for mito)'},
    {'ID': 'gnomad_HomAlt', 'Number': '1', 'Type': 'Integer', 'Description': 'gnomAD HomAlt (unpopulated for mito)'},
    {'ID': 'gene_id', 'Number': '1', 'Type': 'String', 'Description': 'Green gene this row is labelled against'},
    {'ID': 'csq', 'Number': '.', 'Type': 'String', 'Description': 'Talos-formatted transcript consequences'},
]


def write_empty_vcf(vcf_path: str, output_path: str):
    """Write a zero-variant copy of the input VCF - used when there is no mito analysis to do."""
    reader = VCF(vcf_path)
    writer = Writer(output_path, reader, mode='wz')
    writer.close()
    reader.close()
    logger.info(f'Wrote zero-variant VCF to {output_path}')


def cli_main():
    parser = ArgumentParser(description='Labels a BCSQ-annotated mito VCF with fresh ClinVar decisions')
    parser.add_argument('--input', help='Path to the annotated mito VCF', required=True)
    parser.add_argument('--output', help='output VCF path, written bgzipped', required=True)
    parser.add_argument('--panelapp', help='Path to a PanelApp data file for the cohort', required=True)
    parser.add_argument('--pedigree', help='Path to the pedigree file for the cohort', required=True)
    parser.add_argument('--clinvar', help='Path to a ClinvArbitration decisions TSV', required=True)
    args = parser.parse_args()

    main(
        vcf_path=args.input,
        output_path=args.output,
        panelapp=args.panelapp,
        pedigree=args.pedigree,
        clinvar_path=args.clinvar,
    )


def main(
    vcf_path: str,
    output_path: str,
    panelapp: str,
    pedigree: str,
    clinvar_path: str,
):
    """
    Stream the annotated mito VCF, labelling ClinVar P/LP variants in green mito genes.

    Args:
        vcf_path (str): path to the annotated mito VCF
        output_path (str): path to write the labelled VCF to
        panelapp (str): path to the panelapp data file
        pedigree: path to the pedigree file
        clinvar_path (str): path to the ClinvArbitration decisions TSV
    """

    panel_data = read_json_from_path(panelapp, return_model=PanelApp)

    # gets a lookup of all relevant genes from the PanelApp dump
    symbol_to_ensg = {gene.symbol: key for key, gene in panel_data.genes.items() if gene.chrom.startswith('M')}

    # there's no mito analysis to do if there are no mito genes on the panel
    if not symbol_to_ensg:
        logger.info('No mito genes on the panel, no analysis to do')
        write_empty_vcf(vcf_path, output_path)
        sys.exit(0)

    green_ensg_ids = set(symbol_to_ensg.values())

    pedigree_data = parse_pedigree(pedigree)

    reader = VCF(vcf_path)

    # subset to currently considered samples
    if not subset_reader_to_pedigree(reader, pedigree_data):
        logger.error('No samples shared between pedigree and VCF')
        reader.close()
        write_empty_vcf(vcf_path, output_path)
        sys.exit(0)

    # load the fresh ClinVar decisions, restricted to mito contigs - a tiny lookup
    clinvar_decisions = read_clinvar_decisions_tsv(clinvar_path, contigs=MITO_CONTIGS)

    # pull and split the CSQ header line
    csq_fields = split_csq_header(reader)

    for header_entry in OUTPUT_INFO_HEADERS:
        reader.add_info_to_header(header_entry)

    writer = Writer(output_path, reader, mode='wz')

    kept = 0
    for variant in reader:
        # remove filter-failed variants
        if not variant_is_pass(variant):
            continue

        # only ClinVar Pathogenic with 1+ stars survives the mito process
        decision = clinvar_decisions.get(
            (normalise_chrom(variant.CHROM), variant.POS, variant.REF, variant.ALT[0]),
        )
        if decision is None or decision['clinical_significance'] != PATHOGENIC or decision['gold_stars'] < 1:
            continue

        bcsq = variant.INFO.get('BCSQ')
        if not bcsq:
            continue

        consequences = parse_bcsq_entries(bcsq, csq_fields)
        for csq in consequences:
            # work backwards, use panelapp to get ENSGs. BCFtools isn't providing this (no MT MANE data)
            csq['gene_id'] = symbol_to_ensg.get(csq['gene'], csq['gene'])

        green_gene_ids = {csq['gene_id'] for csq in consequences} & green_ensg_ids
        if not green_gene_ids:
            continue

        del variant.INFO['BCSQ']
        variant.INFO['clinvar_significance'] = decision['clinical_significance']
        variant.INFO['clinvar_stars'] = decision['gold_stars']
        variant.INFO['clinvar_allele'] = decision['allele_id']
        variant.INFO['clinvar_talos'] = 1
        variant.INFO['categorybooleanclinvarplp'] = 1
        # todo generate gnomAD Mito annotations from some source
        variant.INFO['gnomad_AC'] = 0
        variant.INFO['gnomad_AF'] = 0.0
        variant.INFO['gnomad_AC_XY'] = 0
        variant.INFO['gnomad_HomAlt'] = 0
        variant.INFO['csq'] = consequences_to_csq_string(consequences)

        # one output row per green gene this variant sits in
        for gene_id in sorted(green_gene_ids):
            variant.INFO['gene_id'] = gene_id
            writer.write_record(variant)
            kept += 1

    writer.close()
    reader.close()
    logger.info(f'Wrote {kept} labelled mito rows to {output_path}')


if __name__ == '__main__':
    cli_main()
