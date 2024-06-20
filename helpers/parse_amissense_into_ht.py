#! /usr/bin/env python3

"""
Script for parsing the compressed tsv file of AlphaMissense results into a Hail Table
This Hail Table can be used for annotation of an old VEP'd VCF

This is being included as annotation using AlphaMissense is important for Talos
We've removed attempts to find consensus using polyphen, SIFT, Revel, CADD, etc.
finding instead that there was near-perfect overlap between the results when using
multiple in Silico tools vs. just AlphaMissense

The AlphaMissense data is publicly available, and this is just a parser for that TSV

Assumption: the first lines of the AM data are still
# Copyright 2023 DeepMind Technologies Limited
#
# Licensed under CC BY-NC-SA 4.0 license

This script makes no attempt to modify or subvert that license, merely rearrange the data

Process:
1. read through the compressed data, skip non-pathogenic entries
2. write the pathogenic entries back out to a reduced schema
3. parse that data as a Hail Table using specific variant types
4. write the Hail Table
"""

import gzip
import os
import string
from argparse import ArgumentParser
from random import choices

import hail as hl


def filter_for_pathogenic_am(input_file: str, intermediate_file: str):
    """
    read the tsv file, skim for pathogenic entries, then write out to a new file

    Args:
        input_file ():
        intermediate_file ():
    """

    headers = ['chrom', 'pos', 'ref', 'alt', 'transcript', 'am_pathogenicity', 'am_class']
    with gzip.open(input_file, 'rt') as read_handle:
        with gzip.open(intermediate_file, 'wt') as write_handle:
            write_handle.write('\t'.join(headers) + '\n')
            for line in read_handle:
                # skip over the headers
                if line.startswith('#'):
                    continue

                # skip over everything except pathogenic
                if 'pathogenic' not in line:
                    continue

                content = line.rstrip().split()
                new_content = [
                    content[0],  # chrom
                    content[1],  # position
                    content[2],  # ref
                    content[3],  # alt
                    content[6].split('.')[0],  # transcript, with version decimal removed
                    content[8],  # am_pathogenicity
                    content[9],  # am_class
                ]

                write_handle.write('\t'.join(new_content) + '\n')


def hail_table_from_tsv(tsv_file: str, new_ht: str):
    """
    take the TSV file previously created and ingest it as a Hail Table
    Args:
        tsv_file ():
        new_ht ():
    """
    hl.init()
    hl.default_reference('GRCh38')
    # import as a hail table, force=True as this isn't Block-Zipped so all read on one core
    # not too bad, it's a small table.
    # We also provide some data types for non-string columns
    ht = hl.import_table(tsv_file, types={'am_pathogenicity': hl.tfloat64, 'pos': hl.tint32}, force=True)

    # combine the two alleles into a single list
    ht = ht.transmute(locus=hl.locus(contig=ht.chrom, pos=ht.pos), alleles=[ht.ref, ht.alt])
    ht = ht.key_by('locus', 'alleles')
    ht.write(new_ht)


if __name__ == '__main__':
    # two positional
    parser = ArgumentParser()
    parser.add_argument('am_tsv', help='path to the AM tsv.gz file')
    parser.add_argument('ht_out', help='path to write a new Hail Table')
    args, unknown = parser.parse_known_args()

    if unknown:
        raise ValueError(unknown)

    # generate a random file name so that we don't overwrite anything consistently
    random_intermediate_file: str = ''.join(choices(string.ascii_uppercase + string.digits, k=6)) + '.tsv.gz'

    # generate a new tsv of just pathogenic entries
    filter_for_pathogenic_am(args.am_tsv, random_intermediate_file)

    # now ingest as HT and re-jig some fields
    hail_table_from_tsv(random_intermediate_file, args.ht_out)

    # if that succeeded, delete the intermediate file
    os.remove(random_intermediate_file)
