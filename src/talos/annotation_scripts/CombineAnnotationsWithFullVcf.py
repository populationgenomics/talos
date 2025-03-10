#!/usr/bin/env python3

"""
Takes the VCF joint-call, and the sites-only annotations reformatted as a HailTable, and hops the annotations from one
to the other.
"""

import json
import logging
from argparse import ArgumentParser
from collections import defaultdict

import hail as hl


MISSING_STRING = hl.str('')


def extract_and_split_csq_string(vcf_path: str) -> list[str]:
    """
    Extract the BCSQ header from the VCF and split it into a list of strings

    Args:
        vcf_path (str): path to the local VCF

    Returns:
        list of strings
    """

    # get the headers from the VCF
    all_headers = hl.get_vcf_metadata(vcf_path)

    # get the '|'-delimited String of all header names
    csq_whole_string = all_headers['info']['BCSQ']['Description'].split('Format: ')[-1]

    # split it all on pipes, return the list
    return csq_whole_string.lower().split('|')


def csq_strings_into_hail_structs(csq_strings: list[str], mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Take the list of BCSQ strings, split the CSQ annotation and re-organise as a hl struct

    Args:
        csq_strings (list[str]): a list of strings, each representing a CSQ entry
        mt (hl.MatrixTable): the MatrixTable to annotate

    Returns:
        a MatrixTable with the BCSQ annotations re-arranged
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
            hl.len(x) == 4,  # noqa: PLR2004
            x.extend([MISSING_STRING, MISSING_STRING, MISSING_STRING]),
            hl.if_else(
                # 5 values... add 2 missing Strings
                hl.len(x) == 5,  # noqa: PLR2004
                x.extend([MISSING_STRING, MISSING_STRING]),
                hl.if_else(
                    hl.len(x) == 6,  # noqa: PLR2004
                    x.extend([MISSING_STRING]),
                    x,
                ),
            ),
        ),
    )

    # transform the CSQ string arrays into structs using the header names
    # Consequence | gene | transcript | biotype | strand | amino_acid_change | dna_change
    mt = mt.annotate_rows(
        transcript_consequences=split_csqs.map(
            lambda x: hl.struct(
                **{csq_strings[n]: x[n] for n in range(len(csq_strings))},
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
            ),
            mt.transcript_consequences,
        ),
    )


def annotate_gene_ids(mt: hl.MatrixTable, bed_file: str) -> hl.MatrixTable:
    """
    The BED file contains the gene IDs, but not all is applied by BCFtool csq
    This method will add the gene IDs to the MatrixTable

    Args:
        mt ():
        bed_file (str): path to a bed file containing gene IDs
    """

    # indexed on contig, then gene symbol: ID
    id_dict: dict[str, dict[str, str]] = defaultdict(dict)

    with open(bed_file) as handle:
        for line in handle:
            # skip over headers and dividing lines
            if line.startswith('#'):
                continue

            chrom, _start, _end, details = line.rstrip().split('\t')
            ensg, symbol = details.split(';')
            id_dict[chrom][symbol] = ensg

    id_hl_dict = hl.literal(id_dict)

    # take the ENSG value from the dict for the contig (correctly matches PAR region genes)
    # deafult to the gene symbol (which can be the ENSG, depending on transcript consequence)
    return mt.annotate_rows(
        transcript_consequences=hl.map(
            lambda x: x.annotate(
                gene_id=id_hl_dict[mt.locus.contig].get(x.gene, x.gene),
            ),
            mt.transcript_consequences,
        ),
    )


def insert_am_annotations(mt: hl.MatrixTable, am_table: str | None = None) -> hl.MatrixTable:
    """
    Load up a Hail Table of AlphaMissense annotations, and annotate this data unless the AM annotations already exist

    Args:
        mt ():
        am_table (str | None):
    """

    if am_table is None:
        logging.error('AM annotations table is not present, and AM annotations are not in the VCF. Please create')

    logging.info(f'Reading AM annotations from {am_table} and applying to MT')

    # read in the hail table containing alpha missense annotations
    am_ht = hl.read_table(am_table)

    # gross - this needs a conditional application based on the specific transcript
    return mt.annotate_rows(
        transcript_consequences=hl.map(
            lambda x: x.annotate(
                am_class=hl.if_else(
                    x.transcript == am_ht[mt.row_key].transcript,
                    am_ht[mt.row_key].am_class,
                    MISSING_STRING,
                ),
                am_pathogenicity=hl.if_else(
                    x.transcript == am_ht[mt.row_key].transcript,
                    am_ht[mt.row_key].am_pathogenicity,
                    hl.float64(0),
                ),
            ),
            mt.transcript_consequences,
        ),
    )


def apply_mane_annotations(mt: hl.MatrixTable, mane_path: str | None = None) -> hl.MatrixTable:
    """
    Apply MANE annotations to the VCF

    Args:
        mt ():
        mane_path (str | None): path to a Hail Table containing MANE annotations

    Returns:
        The same MatrixTable but with additional annotations
    """

    if mane_path is None:
        logging.info('No MANE table found, skipping annotation')
        return mt.annotate_rows(
            transcript_consequences=hl.map(
                lambda x: x.annotate(
                    mane=MISSING_STRING,
                    ensp=MISSING_STRING,
                    nm_id=MISSING_STRING,
                ),
                mt.transcript_consequences,
            ),
        )

    # read in the mane table
    with open(mane_path) as handle:
        mane_dict = json.load(handle)

    # convert the dict into a Hail Dict
    hl_mane_dict = hl.dict(mane_dict)

    key_set = hl_mane_dict.key_set()

    # annotate the variant with MANE data, recovering Mane status, NM ID, and ENSP ID
    return mt.annotate_rows(
        transcript_consequences=hl.map(
            lambda x: x.annotate(
                mane=hl.if_else(key_set.contains(x.transcript), hl_mane_dict[x.transcript]['mane'], MISSING_STRING),
                ensp=hl.if_else(key_set.contains(x.transcript), hl_mane_dict[x.transcript]['ensp'], MISSING_STRING),
                nm_id=hl.if_else(key_set.contains(x.transcript), hl_mane_dict[x.transcript]['nm_id'], MISSING_STRING),
            ),
            mt.transcript_consequences,
        ),
    )


def cli_main():
    """
    take an input VCF and an output MT path
    also supply the alpha_missense table created by parse_amissense_into_ht.py
    """

    parser = ArgumentParser(description='Takes a BCSQ annotated VCF and makes it a MT')
    parser.add_argument('--input', help='Path to the annotated VCF')
    parser.add_argument('--output', help='output MatrixTable path')
    parser.add_argument('--gene_bed', help='BED file containing gene mapping')
    parser.add_argument('--am', help='Hail Table containing AlphaMissense annotations', default=None)
    parser.add_argument('--mane', help='Hail Table containing MANE annotations', default=None)
    args, unknown = parser.parse_known_args()
    if unknown:
        raise ValueError(f'Whats the deal with {unknown}?')

    main(vcf_path=args.input, output_path=args.output, gene_bed=args.gene_bed, alpha_m=args.am, mane=args.mane)


def main(vcf_path: str, output_path: str, gene_bed: str, alpha_m: str | None = None, mane: str | None = None):
    """
    Takes a VEP-annotated VCF, reorganises into a Talos-compatible MatrixTable
    If supplied, will annotate at runtime with AlphaMissense annotations

    Args:
        vcf_path ():
        output_path ():
        gene_bed (str):
        alpha_m ():
        mane (str | None): path to a MANE Hail Table for enhanced annotation
    """

    hl.default_reference('GRCh38')

    # pull and split the CSQ header line
    csq_fields = extract_and_split_csq_string(vcf_path=vcf_path)

    # read the VCF into a MatrixTable
    mt = hl.import_vcf(vcf_path, array_elements_required=False, force_bgz=True)

    # checkpoint it locally to make everything faster
    mt = mt.checkpoint('checkpoint.mt', overwrite=True, _read_if_exists=True)

    # re-shuffle the BCSQ elements
    mt = csq_strings_into_hail_structs(csq_fields, mt)

    # add ENSG IDs where possible
    mt = annotate_gene_ids(mt, bed_file=gene_bed)

    # get a hold of the geneIds - use some aggregation
    mt = mt.annotate_rows(gene_ids=hl.set(mt.transcript_consequences.map(lambda c: c.gene_id)))

    # add AlphaMissense scores
    mt = insert_am_annotations(mt, am_table=alpha_m)

    # drop the BCSQ field
    mt = mt.annotate_rows(info=mt.info.drop('BCSQ'))

    # checkpoint again before rummaging around in MANE table
    mt = mt.checkpoint('checkpoint_pre_mane.mt', overwrite=True, _read_if_exists=True)

    mt = apply_mane_annotations(mt, mane_path=mane)

    mt.describe()

    mt.write(output_path, overwrite=True)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    cli_main()
