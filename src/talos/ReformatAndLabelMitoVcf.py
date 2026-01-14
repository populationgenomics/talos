#!/usr/bin/env python3

"""
This is an adapter process to take a sites-only VCF annotated with gnomAD frequencies and BCFtools csq consequences, and
re-arrange it into a HailTable for use with the Talos pipeline.

This process combines the AF/CSQs already applied with the MANE transcript/protein names, and AlphaMissense annotations
"""

import sys
from argparse import ArgumentParser

from loguru import logger
from mendelbrot.pedigree_parser import PedigreeParser

import hail as hl

from cpg_utils.hail_batch import init_batch

from talos.annotation_scripts.ReformatAnnotatedVcfIntoHailTable import (
    MISSING_FLOAT,
    MISSING_INT,
    MISSING_STRING,
    extract_and_split_csq_string,
)
from talos.config import config_retrieve
from talos.models import PanelApp
from talos.RunHailFiltering import (
    annotate_clinvarbitration,
    csq_struct_to_string,
    filter_matrix_by_ac,
    subselect_mt_to_pedigree,
)
from talos.utils import read_json_from_path


def csq_strings_into_hail_structs(csq_strings: list[str], mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Take the list of BCSQ strings, split the CSQ annotation and re-organise as a hl struct

    Args:
        csq_strings (list[str]): a list of strings, each representing a CSQ entry
        mt (hl.Table): the Table to annotate

    Returns:
        a Table with the BCSQ annotations re-arranged
    """

    # get the BCSQ contents as a list of lists of strings, per variant
    split_csqs = mt.info.BCSQ.map(lambda csq_entry: csq_entry.split('\|'))  # noqa: W605

    # this looks pretty hideous, bear with me
    # if BCFtools csq doesn't have a consequence annotation, it will truncate the pipe-delimited string
    # this is fine sometimes, but not when we're building a schema here
    # when we find truncated BCSQ strings, we need to add dummy values to the end of the array
    split_csqs = split_csqs.map(
        lambda x: hl.if_else(
            # if there were only 4 values, add 3 missing Strings
            hl.len(x) == 4,
            x.extend([MISSING_STRING, MISSING_STRING, MISSING_STRING]),
            hl.if_else(
                # 5 values... add 2 missing Strings
                hl.len(x) == 5,
                x.extend([MISSING_STRING, MISSING_STRING]),
                hl.if_else(
                    hl.len(x) == 6,
                    x.extend([MISSING_STRING]),
                    hl.if_else(
                        hl.len(x) == 3,
                        x.extend([MISSING_STRING, MISSING_STRING, MISSING_STRING, MISSING_STRING]),
                        x,
                    ),
                ),
            ),
        ),
    )

    # transform the CSQ string arrays into structs using the header names
    # Consequence | gene | transcript | biotype | strand | amino_acid_change | dna_change
    mt = mt.annotate_rows(
        transcript_consequences=split_csqs.map(
            lambda x: hl.struct(
                **{csq_strings[n]: x[n] for n in range(len(csq_strings)) if csq_strings[n] != 'strand'},
            ),
        ),
    )

    return mt.annotate_rows(
        # amino_acid_change can be absent, or in the form of "123P" or "123P-124F"
        # we use this number when matching to the codons of missense variants, to find codon of the reference pos.
        transcript_consequences=hl.map(
            lambda x: x.annotate(
                codon=hl.if_else(
                    x.amino_acid_change == MISSING_STRING,
                    hl.missing(hl.tint32),
                    hl.if_else(
                        x.amino_acid_change.matches('^([0-9]+).*$'),
                        hl.int32(x.amino_acid_change.replace('^([0-9]+).+', '$1')),
                        hl.missing(hl.tint32),
                    ),
                ),
                am_class=MISSING_STRING,
                am_pathogenicity=MISSING_FLOAT,
                # MANE transcript selection has not been applied to mitochondrial data
                mane_status=MISSING_STRING,
                ensp=MISSING_STRING,
                mane_id=MISSING_STRING,
            ),
            mt.transcript_consequences,
        ),
    )


def cli_main():
    """
    take an input VCF and an output MT path, apply the correct reformatting/schema to the annotations
    """

    parser = ArgumentParser(description='Takes a BCSQ annotated VCF and makes it a mt')
    parser.add_argument('--input', help='Path to the annotated sites-only VCF', required=True)
    parser.add_argument('--output', help='output Table path, must have a ".mt" extension', required=True)
    parser.add_argument('--panelapp', help='Path to a PanelApp data file for the cohort', required=True)
    parser.add_argument('--pedigree', help='Path to the pedigree file for the cohort', required=True)
    parser.add_argument('--clinvar', help='Path to a ClinvArbitration decisions HT', required=True)
    parser.add_argument('--batch', help='flag to use the batch hail backend', action='store_true')
    args = parser.parse_args()

    main(
        vcf_path=args.input,
        output_path=args.output,
        panelapp=args.panelapp,
        pedigree=args.pedigree,
        clinvar_path=args.clinvar,
        batch=args.batch,
    )


def main(
    vcf_path: str,
    output_path: str,
    panelapp: str,
    pedigree: str,
    clinvar_path: str,
    batch: bool,
):
    """
    Takes a BCFtools-annotated Mito VCF, reorganises into a Talos-compatible MatrixTable

    Args:
        vcf_path (str): path to the annotated sites-only VCF
        output_path (str): path to write the resulting Hail Table to
        panelapp (str): path to the panelapp data file
        pedigree: path to the pedigree file
        clinvar_path (str): path to the clinvar data file
        batch (bool): if we should use the Hail Batch backend
    """

    panel_data = read_json_from_path(panelapp, return_model=PanelApp)

    # gets a lookup of all relevant genes from the PanelApp dump
    symbol_to_ensg = {gene.symbol: key for key, gene in panel_data.genes.items() if gene.chrom.startswith('M')}

    if batch:
        logger.info('Using Hail Batch backend')
        init_batch()
    else:
        logger.info('Using Hail Local backend')
        hl.context.init_spark(master='local[*]', default_reference='GRCh38', quiet=True)

    # pull and split the CSQ header line
    csq_fields = extract_and_split_csq_string(vcf_path=vcf_path)

    # read the VCF into a MatrixTable
    mt = hl.import_vcf(vcf_path, array_elements_required=False, force_bgz=True)

    # there's no mito analysis to do if there are no mito genes on the panel
    # this will be resolved by an enhanced mendeliome at some point, but for now we need to quit out
    # drop all the rows, and write a zero-variant VCF. It's a messy compromise, I'll fix it properly later
    if not symbol_to_ensg:
        mt = mt.filter_rows(hl.is_defined(mt.filters), keep=False)
        hl.export_vcf(mt, output_path, tabix=True)
        sys.exit(0)

    symbol_to_ensg_hl = hl.literal(symbol_to_ensg)
    green_ensg_ids = hl.literal(set(symbol_to_ensg.values()))

    # remove filtered variants
    mt = mt.filter_rows(hl.is_missing(mt.filters) | (mt.filters.length() == 0))

    # parse the pedigree into an object
    pedigree_data = PedigreeParser(pedigree)

    # reduce cohort to singletons, if the config says so
    if config_retrieve('singletons', False):
        logger.info('Reducing pedigree to affected singletons only')
        pedigree_data.set_participants(pedigree_data.as_singletons())
        pedigree_data.set_participants(pedigree_data.get_affected_members())

    # subset to currently considered samples
    try:
        mt = subselect_mt_to_pedigree(mt, ped_samples=pedigree_data.get_all_sample_ids())
    except ValueError as ve:
        logger.error(str(ve))
        mt = mt.filter_rows(hl.is_defined(mt.filters), keep=False)
        hl.export_vcf(mt, output_path, tabix=True)
        sys.exit(0)

    mt = annotate_clinvarbitration(mt=mt, clinvar=clinvar_path)

    mt = filter_matrix_by_ac(mt)

    # re-shuffle the BCSQ elements
    mt = csq_strings_into_hail_structs(csq_fields, mt)

    mt = mt.annotate_rows(
        # work backwards, use panelapp to get ENSGs. BCFtools isn't providing this, we get it from MANE (No MT MANE)
        transcript_consequences=hl.map(
            lambda x: x.annotate(
                gene_id=symbol_to_ensg_hl.get(x.gene, x.gene),
            ),
            mt.transcript_consequences,
        ),
        # todo generate gnomAD Mito annotations from some source
        gnomad=hl.struct(
            gnomad_AC=MISSING_INT,
            gnomad_AF=MISSING_FLOAT,
            gnomad_AC_XY=MISSING_INT,
            gnomad_HomAlt=MISSING_INT,
        ),
    )

    # get a hold of the geneIds - use some aggregation
    mt = mt.annotate_rows(gene_ids=hl.set(mt.transcript_consequences.map(lambda c: c.gene_id)))

    # drop the BCSQ field
    mt = mt.annotate_rows(info=mt.info.drop('BCSQ'))

    mt = mt.explode_rows(mt.gene_ids)

    # hard filter remaining rows to PanelApp green genes
    mt = mt.filter_rows(green_ensg_ids.contains(mt.gene_ids))

    # filter everything else out
    mt = mt.filter_rows(mt.info.categorybooleanclinvarplp == 1)

    # obtain the massive CSQ string using method stolen from the Broad's Gnomad library
    # also take the single gene_id (from the exploded attribute)
    # and retrieves the gnomAD annotations
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            **mt.gnomad,
            csq=csq_struct_to_string(mt.transcript_consequences),
            gene_id=mt.gene_ids,
        ),
    )

    # now write that badboi
    hl.export_vcf(mt, output_path, tabix=True)


if __name__ == '__main__':
    cli_main()
