#!/usr/bin/env python3

"""
Rename SVAFotate's population frequency annotations to the field names Talos expects.

SVAFotate writes `Max_AF`, and with `-a best` also `Best_gnomAD_ID` and `gnomAD_Count`.
`talos.run_hail_filtering_sv` reads `{prefix}_sv_AF` and `{prefix}_sv_SVID`. The prefix is read from the
same config key the filtering script uses, so the two cannot drift apart.

The original SVAFotate fields are left in place as provenance.

Two details of SVAFotate's output are load-bearing here:

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

# the fields run_hail_filtering_sv reads
TALOS_AF = f'{GNOMAD_POP}_sv_AF'
TALOS_SVID = f'{GNOMAD_POP}_sv_SVID'


def rename_annotations(input_path: str, output_path: str):
    """
    Copy SVAFotate's AF and SVID annotations into the field names Talos reads

    Args:
        input_path (str): VCF annotated by SVAFotate
        output_path (str): where to write the renamed VCF
    """

    vcf = VCF(input_path)

    vcf.add_info_to_header(
        {
            'ID': TALOS_AF,
            'Number': '1',
            'Type': 'Float',
            'Description': f'{GNOMAD_POP} SV allele frequency, from SVAFotate {SVAFOTATE_AF}',
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

    writer = Writer(output_path, vcf)

    variants = 0
    with_af = 0
    with_svid = 0

    for variant in vcf:
        variants += 1

        max_af = variant.INFO.get(SVAFOTATE_AF)
        if max_af is not None:
            variant.INFO[TALOS_AF] = max_af
            if max_af > 0:
                with_af += 1

        # only carry the SVID across if a match actually cleared the overlap threshold
        if variant.INFO.get(SVAFOTATE_COUNT):
            best_id = variant.INFO.get(SVAFOTATE_SVID)
            if best_id is not None:
                variant.INFO[TALOS_SVID] = best_id
                with_svid += 1

        writer.write_record(variant)

    writer.close()
    vcf.close()

    logger.info(f'Renamed annotations on {variants} variants')
    logger.info(f'{with_af} carried a non-zero {TALOS_AF}, {with_svid} carried a {TALOS_SVID}')

    if variants and not with_af:
        logger.warning(
            f'No variant received a non-zero {TALOS_AF}. Every variant will look rare to downstream '
            'filtering. Check that the SVAFotate reference BED uses Ensembl-style contig names (1, not chr1) '
            'and that the requested source was present in it.',
        )


def cli_main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, help='VCF annotated by SVAFotate')
    parser.add_argument('--output', required=True, help='Path to write the renamed VCF')
    args = parser.parse_args()

    rename_annotations(input_path=args.input, output_path=args.output)


if __name__ == '__main__':
    cli_main()
