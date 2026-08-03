#!/usr/bin/env python3

"""
Rename SVAFotate's population-frequency INFO fields to the names talos.run_hail_filtering_sv reads.

SVAFotate writes `Max_AF`, and with `-a best` also `Best_gnomAD_ID` and `gnomAD_Count`.
`talos.run_hail_filtering_sv` reads `{prefix}_sv_AF` and `{prefix}_sv_SVID`, where `{prefix}` is
`RunHailFilteringSv.gnomad_population`. Reading that prefix from the same config key both sides use keeps the
written and expected field names from drifting apart.

The original SVAFotate fields are left in place as provenance.

Three details of SVAFotate's output matter here:

- use `Max_AF`, not `Best_gnomAD_AF`. `Max_AF` is the maximum frequency across all qualifying matches, which
  is the conservative choice for rare-disease filtering - under-filtering a common variant is a worse outcome
  than losing a rare one.
- `Best_gnomAD_ID` reports the highest-ranked *candidate* overlap whether or not it cleared the reciprocal
  overlap threshold passed to `svafotate --minf`. Copying it unconditionally would attach a gnomAD SVID to
  nearly every variant in the callset, so the SVID is only carried across when `gnomAD_Count` shows a
  qualifying match was actually found.
- SVAFotate matches on coordinate overlap alone, and strips the `chr` prefix from the query VCF but not from
  the reference BED. A contig naming disagreement therefore yields `Max_AF=0` for every variant with no
  error raised, which looks exactly like a callset of rare variants. If nothing at all matched, this script
  warns loudly rather than passing the result on silently.
"""

from argparse import ArgumentParser

from cyvcf2 import VCF, Writer
from loguru import logger

from talos.config import config_retrieve

GNOMAD_POP = config_retrieve(['RunHailFilteringSv', 'gnomad_population'], 'gnomad_v4.1')

# the fields SVAFotate writes, and the source whose count gates the SVID
SVAFOTATE_AF = 'Max_AF'
SVAFOTATE_SVID = 'Best_gnomAD_ID'
SVAFOTATE_COUNT = 'gnomAD_Count'

# the names run_hail_filtering_sv reads, derived from the configured population prefix
TALOS_AF = f'{GNOMAD_POP}_sv_AF'
TALOS_SVID = f'{GNOMAD_POP}_sv_SVID'


def cli_main():
    """
    main method wrapper for console script execution
    """
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, help='SVAFotate-annotated VCF')
    parser.add_argument('--output', required=True, help='Where to write the renamed VCF')
    args = parser.parse_args()
    main(vcf_in=args.input, vcf_out=args.output)


def main(vcf_in: str, vcf_out: str):
    """
    Copy SVAFotate's Max_AF/Best_gnomAD_ID into the field names Talos reads, and write the VCF back out.

    Args:
        vcf_in (str): path to the SVAFotate-annotated VCF
        vcf_out (str): path to write the renamed VCF to
    """

    vcf = VCF(vcf_in)

    vcf.add_info_to_header(
        {
            'ID': TALOS_AF,
            'Number': '1',
            'Type': 'Float',
            'Description': f'{GNOMAD_POP} SV allele frequency, copied from SVAFotate {SVAFOTATE_AF}',
        },
    )
    vcf.add_info_to_header(
        {
            'ID': TALOS_SVID,
            'Number': '1',
            'Type': 'String',
            'Description': f'{GNOMAD_POP} SV ID of the best qualifying match, from SVAFotate {SVAFOTATE_SVID}',
        },
    )

    writer = Writer(vcf_out, vcf)

    variants = 0
    with_af = 0
    with_svid = 0

    for variant in vcf:
        variants += 1

        # a variant SVAFotate left without Max_AF had no candidate at all, which is frequency zero
        max_af = variant.INFO.get(SVAFOTATE_AF, 0.0)
        variant.INFO[TALOS_AF] = max_af
        if max_af > 0:
            with_af += 1

        # only carry the SVID across if a match actually cleared the overlap threshold
        if variant.INFO.get(SVAFOTATE_COUNT, 0) > 0 and (svid := variant.INFO.get(SVAFOTATE_SVID)):
            variant.INFO[TALOS_SVID] = svid
            with_svid += 1

        writer.write_record(variant)

    writer.close()
    vcf.close()

    logger.info(f'Renamed annotations on {variants} variants')
    logger.info(f'{with_af} carried a non-zero {TALOS_AF}, {with_svid} carried a {TALOS_SVID}')

    if variants and not with_af:
        logger.warning(
            f'No variant carried a non-zero {SVAFOTATE_AF}. Every variant will look rare to downstream '
            'filtering. Check that the SVAFotate reference BED uses Ensembl-style contig names (1, not chr1) '
            'and that the requested source was present in it.',
        )


if __name__ == '__main__':
    cli_main()
