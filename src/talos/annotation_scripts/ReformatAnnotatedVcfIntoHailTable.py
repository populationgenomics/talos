#!/usr/bin/env python3

"""
This is an adapter process to take a sites-only VCF annotated with gnomAD frequencies and BCFtools csq consequences, and
re-arrange it into a HailTable for use with the Talos pipeline.

This process combines the AF/CSQs already applied with the MANE transcript/protein names, and AlphaMissense annotations
"""

import json
from argparse import ArgumentParser
from collections import defaultdict

from loguru import logger

import hail as hl

from cpg_utils.hail_batch import init_batch

MISSING_FLOAT = hl.float64(0)
MISSING_INT = hl.int32(0)
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


def csq_strings_into_hail_structs(csq_strings: list[str], ht: hl.Table | hl.MatrixTable) -> hl.Table | hl.MatrixTable:
    """
    Take the list of BCSQ strings, split the CSQ annotation and re-organise as a hl struct

    Args:
        csq_strings (list[str]): a list of strings, each representing a CSQ entry
        ht (hl.Table): the Table to annotate

    Returns:
        a Table with the BCSQ annotations re-arranged
    """

    # get the BCSQ contents as a list of lists of strings, per variant
    split_csqs = ht.info.BCSQ.map(lambda csq_entry: csq_entry.split('\|'))  # noqa: W605

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
                    x,
                ),
            ),
        ),
    )

    # transform the CSQ string arrays into structs using the header names
    # Consequence | gene | transcript | biotype | strand | amino_acid_change | dna_change
    ht = ht.annotate(
        transcript_consequences=split_csqs.map(
            lambda x: hl.struct(
                **{csq_strings[n]: x[n] for n in range(len(csq_strings)) if csq_strings[n] != 'strand'},
            ),
        ),
    )

    return ht.annotate(
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
            ht.transcript_consequences,
        ),
    )


def annotate_gene_ids(ht: hl.Table, bed_file: str) -> hl.Table:
    """
    Using a BED file to generate lookup, annotate each transcript consequence with the Ensembl gene ID(s).
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
    # default to the gene symbol (which can be the ENSG, depending on transcript consequence)
    return ht.annotate(
        transcript_consequences=hl.map(
            lambda x: x.annotate(
                gene_id=id_hl_dict[ht.locus.contig].get(x.gene, x.gene),
            ),
            ht.transcript_consequences,
        ),
    )


def insert_am_annotations(ht: hl.Table, am_table: str) -> hl.Table:
    """
    Load up a Hail Table of AlphaMissense annotations, and annotate this data unless the AM annotations already exist.
    """

    logger.info(f'Reading AM annotations from {am_table} and applying to MT')

    # read in the hail table containing alpha missense annotations
    am_ht = hl.read_table(am_table)

    # AM consequence matching needs conditional application based on the specific transcript match
    return ht.annotate(
        transcript_consequences=hl.map(
            lambda x: x.annotate(
                am_class=hl.if_else(
                    x.transcript == am_ht[ht.key].transcript,
                    am_ht[ht.key].am_class,
                    MISSING_STRING,
                ),
                am_pathogenicity=hl.if_else(
                    x.transcript == am_ht[ht.key].transcript,
                    am_ht[ht.key].am_pathogenicity,
                    MISSING_FLOAT,
                ),
            ),
            ht.transcript_consequences,
        ),
    )


def apply_mane_annotations(ht: hl.Table, mane_path: str | None = None) -> hl.Table:
    """
    Apply MANE annotations to the VCF.

    If a MANE json file is provided, it will annotate the transcript consequences with the MANE status/ID
    """

    if mane_path is None:
        logger.info('No MANE table found, skipping annotation - dummy values will be entered instead')
        return ht

    # read in the mane table
    with open(mane_path) as handle:
        mane_dict = json.load(handle)

    # convert the dict into a Hail Dict
    hl_mane_dict = hl.dict(mane_dict)

    key_set = hl_mane_dict.key_set()

    # annotate the variant with MANE data, recovering Mane status, NM ID, and ENSP ID
    return ht.annotate(
        transcript_consequences=hl.map(
            lambda x: x.annotate(
                mane_status=hl.if_else(
                    key_set.contains(x.transcript),
                    hl_mane_dict[x.transcript]['mane_status'],
                    MISSING_STRING,
                ),
                ensp=hl.if_else(
                    key_set.contains(x.transcript),
                    hl_mane_dict[x.transcript]['ensp'],
                    MISSING_STRING,
                ),
                mane_id=hl.if_else(
                    key_set.contains(x.transcript),
                    hl_mane_dict[x.transcript]['mane_id'],
                    MISSING_STRING,
                ),
            ),
            ht.transcript_consequences,
        ),
    )


def cli_main():
    """
    take an input VCF and an output MT path
    also supply the alpha_missense table created by parse_amissense_into_ht.py
    """

    parser = ArgumentParser(description='Takes a BCSQ annotated VCF and makes it a HT')
    parser.add_argument('--input', help='Path to the annotated sites-only VCF', required=True)
    parser.add_argument('--output', help='output Table path, must have a ".ht" extension', required=True)
    parser.add_argument('--gene_bed', help='BED file containing gene mapping')
    parser.add_argument('--am', help='Hail Table containing AlphaMissense annotations', required=True)
    parser.add_argument('--mane', help='Hail Table containing MANE annotations', default=None)
    parser.add_argument(
        '--checkpoint',
        help='Whether to use a remote checkpoint. This is an implicit trigger for the batch backend',
        default=None,
    )
    args = parser.parse_args()

    main(
        vcf_path=args.input,
        output_path=args.output,
        gene_bed=args.gene_bed,
        alpha_m=args.am,
        mane=args.mane,
        checkpoint=args.checkpoint,
    )


def main(
    vcf_path: str,
    output_path: str,
    gene_bed: str,
    alpha_m: str,
    mane: str | None = None,
    checkpoint: str | None = None,
):
    """
    Takes a BCFtools-annotated VCF, reorganises into a Talos-compatible MatrixTable
    Will annotate at runtime with AlphaMissense annotations

    Args:
        vcf_path (str): path to the annotated sites-only VCF
        output_path (str): path to write the resulting Hail Table to, must
        gene_bed (str): path to a BED file containing gene IDs, derived from the Ensembl GFF3 file
        alpha_m (str): path to the AlphaMissense Hail Table, required
        mane (str | None): path to a MANE Hail Table for enhanced annotation
        checkpoint (str): which hail backend to use. Defaults to
    """

    if checkpoint:
        logger.info(f'Using Hail Batch backend, checkpointing to {checkpoint}')
        init_batch()
    else:
        logger.info('Using Hail Local backend, no checkpoints')
        hl.context.init_spark(master='local[*]', default_reference='GRCh38', quiet=True)

    # pull and split the CSQ header line
    csq_fields = extract_and_split_csq_string(vcf_path=vcf_path)

    # read the VCF into a MatrixTable
    mt = hl.import_vcf(vcf_path, array_elements_required=False, force_bgz=True)

    # checkpoint the rows as a Table locally to make everything downstream faster
    ht = mt.rows().checkpoint(checkpoint or 'checkpoint.ht', overwrite=True, _read_if_exists=True)

    # re-shuffle the BCSQ elements
    ht = csq_strings_into_hail_structs(csq_fields, ht)

    # add ENSG IDs where possible
    ht = annotate_gene_ids(ht, bed_file=gene_bed)

    # get a hold of the geneIds - use some aggregation
    ht = ht.annotate(gene_ids=hl.set(ht.transcript_consequences.map(lambda c: c.gene_id)))

    # add AlphaMissense scores
    ht = insert_am_annotations(ht, am_table=alpha_m)

    # drop the BCSQ field
    ht = ht.annotate(info=ht.info.drop('BCSQ'))

    ht = apply_mane_annotations(ht, mane_path=mane)

    ht.describe()

    ht.write(output_path, overwrite=True)


if __name__ == '__main__':
    cli_main()
