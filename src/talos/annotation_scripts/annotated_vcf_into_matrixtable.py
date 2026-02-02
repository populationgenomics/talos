#!/usr/bin/env python3

"""
Takes the annotated VCF, samples and all, and reads it as a MatrixTable.
This rearranges all the annotations into the format expected downstream.
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


def csq_strings_into_hail_structs(csq_strings: list[str], mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Take the list of BCSQ strings, split the CSQ annotation and re-organise as a hl struct

    Args:
        csq_strings (list[str]): a list of strings, each representing a CSQ entry
        mt (hl.MatrixTable): the Mt to annotate

    Returns:
        original MatrixTable with the BCSQ annotations re-arranged
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
            ),
            mt.transcript_consequences,
        ),
    )


def get_gene_id_dict(bed_file: str) -> hl.DictExpression:
    """
    Use the BED file to generate a lookup of Gene Symbols to Ensembl Gene IDs
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

    return hl.literal(id_dict)


def get_mane_annotations(mane_path: str) -> hl.DictExpression:
    """
    Parse the MANE file, and get a dict expression of annotations.
    Args:
        mane_path (str): path to the MANE file

    Returns:
        the dict expression of MANE annotations
    """
    # read in the mane table
    with open(mane_path) as handle:
        mane_dict = json.load(handle)

    # convert the dict into a Hail Dict
    return hl.dict(mane_dict)


def annotate_all_transcript_consequences(
    mt: hl.MatrixTable,
    mane: hl.DictExpression,
    ensgs: hl.DictExpression,
) -> hl.MatrixTable:
    """
    In a single loop, update the AM annotations, MANE annotations, and ENSG gene IDs
    Args:
        mt (MatrixTable): the MatrixTable to annotate
        mane (hl.DictExpression): the MANE annotations
        ensgs (hl.DictExpression): the ENSG annotations

    Returns:
        Original MatrixTable, with reformatted + extended annotations
    """
    key_set = mane.key_set()

    return mt.annotate_rows(
        transcript_consequences=hl.map(
            lambda x: x.annotate(
                am_class=hl.if_else(
                    x.transcript == mt.info.am_transcript,
                    mt.info.am_class,
                    MISSING_STRING,
                ),
                am_pathogenicity=hl.if_else(
                    x.transcript == mt.info.am_transcript,
                    mt.info.am_score,
                    MISSING_FLOAT,
                ),
                mane_status=hl.if_else(
                    key_set.contains(x.transcript),
                    mane[x.transcript]['mane_status'],
                    MISSING_STRING,
                ),
                ensp=hl.if_else(
                    key_set.contains(x.transcript),
                    mane[x.transcript]['ensp'],
                    MISSING_STRING,
                ),
                mane_id=hl.if_else(
                    key_set.contains(x.transcript),
                    mane[x.transcript]['mane_id'],
                    MISSING_STRING,
                ),
                gene_id=ensgs[mt.locus.contig].get(x.gene, x.gene),
            ),
            mt.transcript_consequences,
        ),
    )


def nest_gnomad_in_struct(mt: hl.MatrixTable) -> hl.MatrixTable:
    """Tucks all gnomAD annotations into a hl.Struct"""
    return mt.annotate_rows(
        gnomad=hl.struct(
            gnomad_AC=mt.info.gnomad_AC_joint,
            gnomad_AF=mt.info.gnomad_AF_joint,
            gnomad_AC_XY=mt.info.gnomad_AC_joint_XY,
            gnomad_HomAlt=mt.info.gnomad_HomAlt_joint,
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
        mane=args.mane,
        checkpoint=args.checkpoint,
    )


def main(
    vcf_path: str,
    output_path: str,
    gene_bed: str,
    mane: str,
    checkpoint: str | None = None,
):
    """
    Takes a BCFtools-annotated VCF, reorganises into a Talos-compatible MatrixTable
    Will annotate at runtime with AlphaMissense annotations

    Args:
        vcf_path (str): path to the annotated sites-only VCF
        output_path (str): path to write the resulting Hail Table to, must
        gene_bed (str): path to a BED file containing gene IDs, derived from the Ensembl GFF3 file
        mane (str): path to a MANE JSON file for enhanced annotation
        checkpoint (str): which hail backend to use. Defaults to
    """

    if checkpoint:
        logger.info(f'Using Hail Batch backend, checkpointing to {checkpoint}')
        init_batch(
            driver_memory='highmem',
            driver_cores=2,
        )
    else:
        logger.info('Using Hail Local backend, will use a local checkpoint.')
        hl.context.init_spark(master='local[*]', default_reference='GRCh38', quiet=True)

    # pull and split the CSQ header line
    csq_fields = extract_and_split_csq_string(vcf_path=vcf_path)

    # read the VCF into a MatrixTable
    mt = hl.import_vcf(vcf_path, array_elements_required=False, force_bgz=True, block_size=20)

    # checkpoint locally to make everything downstream faster
    mt = mt.checkpoint(checkpoint or 'checkpoint.ht', overwrite=True)

    # re-shuffle the BCSQ elements
    mt = csq_strings_into_hail_structs(csq_fields, mt)

    # Convert JSON data sources into a hl.Dict object
    ensg_dict = get_gene_id_dict(bed_file=gene_bed)
    mane_dict = get_mane_annotations(mane_path=mane)

    # in a single loop, update alphamissense annotations, ENSG gene IDs, and MANE status/matched transcripts
    mt = annotate_all_transcript_consequences(mt, mane_dict, ensg_dict)

    # get a hold of the geneIds - use some aggregation
    mt = mt.annotate_rows(gene_ids=hl.set(mt.transcript_consequences.map(lambda c: c.gene_id)))

    # take note of all named gnomad_* fields
    individual_gnomad_fields = [f for f in mt.info if f.startswith('gnomad_')]

    # gather gnomAD annotations into a separate struct
    mt = nest_gnomad_in_struct(mt)

    # drop the BCSQ field, and all individual gnomAD annotations
    mt = mt.annotate_rows(info=mt.info.drop('BCSQ', *individual_gnomad_fields))

    mt.describe()

    mt.write(output_path, overwrite=True)


if __name__ == '__main__':
    cli_main()
