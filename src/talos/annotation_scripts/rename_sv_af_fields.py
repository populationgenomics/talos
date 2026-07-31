"""
Rename SVAFotate's population-frequency INFO fields to the names talos.run_hail_filtering_sv reads.

SVAFotate writes Max_AF and (with `-a best`) Best_gnomAD_ID, but run_hail_filtering_sv reads
{prefix}_sv_AF and {prefix}_sv_SVID, where {prefix} is RunHailFilteringSv.gnomad_population. Reading that
prefix from the same config key both sides use keeps the written and expected field names from drifting apart.

Two non-obvious rules:
- use Max_AF, not Best_gnomAD_AF - Max_AF is the maximum frequency across all qualifying matches, which is the
  conservative choice for rare-disease filtering
- Best_gnomAD_ID is populated even when nothing cleared the overlap threshold, so the SVID is only carried
  through when gnomAD_Count > 0, and emitted null otherwise
"""

from argparse import ArgumentParser

from cyvcf2 import VCF, Writer
from loguru import logger

from talos.config import config_retrieve

GNOMAD_POP = config_retrieve(['RunHailFilteringSv', 'gnomad_population'], 'gnomad_v4.1')

# SVAFotate's source field names
SVAFOTATE_AF = 'Max_AF'
SVAFOTATE_SVID = 'Best_gnomAD_ID'
SVAFOTATE_COUNT = 'gnomAD_Count'

# the names Talos reads, derived from the configured population prefix
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
            'Description': f'gnomAD-SV population allele frequency, copied from SVAFotate {SVAFOTATE_AF}',
        },
    )
    vcf.add_info_to_header(
        {
            'ID': TALOS_SVID,
            'Number': '1',
            'Type': 'String',
            'Description': f'gnomAD-SV matching variant ID, copied from SVAFotate {SVAFOTATE_SVID}',
        },
    )

    writer = Writer(vcf_out, vcf)

    any_af = False
    for variant in vcf:
        max_af = variant.INFO.get(SVAFOTATE_AF)
        if max_af is None:
            max_af = 0.0
        variant.INFO[TALOS_AF] = max_af
        if max_af > 0:
            any_af = True

        # Best_gnomAD_ID reports the highest-ranked candidate overlap regardless of whether it cleared the
        # threshold, so only carry the SVID through when a qualifying match was actually recorded
        count = variant.INFO.get(SVAFOTATE_COUNT) or 0
        svid = variant.INFO.get(SVAFOTATE_SVID)
        if count > 0 and svid is not None:
            variant.INFO[TALOS_SVID] = svid

        writer.write_record(variant)

    writer.close()
    vcf.close()

    if not any_af:
        logger.warning(
            f'No variant carried a non-zero {SVAFOTATE_AF}. If this is a whole callset, the SVAFotate '
            f'reference BED may have been contig-normalised, which silently zeroes every frequency',
        )


if __name__ == '__main__':
    cli_main()
