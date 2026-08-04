#!/usr/bin/env python3

"""
Sort each complex SV's `CPX_INTERVALS` into coordinate order, so GATK SVAnnotate can process it.

`SVAnnotateEngine.getSegmentForNearestTSS` builds the spanning interval of a complex event by folding its
segments left to right with `SimpleInterval.mergeWithContiguous`, which throws unless each adjacent pair *in
list order* touches. That assumes `CPX_INTERVALS` is coordinate-sorted. GATK-SV does not sort it - it writes
the intervals in the structural order of the event - so a `delINVdel` arrives as

    DEL_chr1:1416597-16725244,DEL_chr1:64416883-234776440,INV_chr1:16725244-64416883

where the two flanking deletions are listed adjacently but are separated by the inversion that bridges them.
The first merge sees the gap and the whole run dies with

    GATKException: The two intervals need to be contiguous: chr1:1416597-16725244 chr1:64416883-234776440

This is a GATK defect (confirmed against 4.6.2.0), not a malformed callset, and it is not recoverable
downstream - SVAnnotate aborts on the first offending record, so a single one costs the entire cohort.

Sorting is safe because order carries no meaning to any other part of SVAnnotate. The gene-overlap, promoter
and non-coding annotators each iterate the segments and accumulate into a dictionary, and the INV adjustment
applied to `dupINV`/`INVdup`/`dupINVdup`/`delINVdup`/`dupINVdel` re-appends its segment at the end regardless
of where it started. The fold in `getSegmentForNearestTSS` is the only order-sensitive step, and there the
sorted order is the one that yields the intended outer-breakpoint span. Talos itself never reads
`CPX_INTERVALS`; `run_hail_filtering_sv` takes `PREDICTED_LOF` and the joint-call fields.

`delINVdel` is the only subtype observed to trigger this - the others reduce to a single segment, or to a
pair that is contiguous by construction, before the fold runs. Every record is sorted anyway rather than
special-casing that one subtype, because the sort is a no-op on an already-ordered record and the cost of
missing a subtype is another aborted cohort.

Sorting cannot rescue a complex event whose segments are genuinely disjoint; SVAnnotate will still reject
those. It only removes the failures caused by the ordering assumption.
"""

from argparse import ArgumentParser

from cyvcf2 import VCF, Writer
from loguru import logger

# a CPX_INTERVALS entry is `SVTYPE_contig:start-end`, e.g. `DEL_chr1:1416597-16725244`
CPX_INTERVALS = 'CPX_INTERVALS'


def interval_sort_key(entry: str, contig_order: dict[str, int]) -> tuple[int, int, int]:
    """
    Build the coordinate sort key for one CPX_INTERVALS entry.

    Contigs are ordered by first appearance rather than by name, so that a multi-contig event keeps its
    segments grouped as the caller wrote them and only the intervals within a contig are reordered. Sorting
    contigs lexically would interleave them, which helps nothing - SVAnnotate cannot merge across contigs
    either way.

    Args:
        entry (str): a single `SVTYPE_contig:start-end` entry
        contig_order (dict[str, int]): contig name -> position of its first appearance in this record

    Returns:
        tuple[int, int, int]: (contig rank, start, end)
    """

    contig, start, end = parse_interval(entry)
    return contig_order[contig], start, end


def parse_interval(entry: str) -> tuple[str, int, int]:
    """
    Split one CPX_INTERVALS entry into its contig and coordinates.

    The SV type is separated from the locus by the first underscore only. GATK's own parser splits on every
    underscore and takes the second field, which mangles contigs such as `chr1_KI270706v1_random`; that is a
    separate GATK bug, and reproducing it here would corrupt the sort rather than match it.

    Args:
        entry (str): a single `SVTYPE_contig:start-end` entry

    Returns:
        tuple[str, int, int]: contig, start, end

    Raises:
        ValueError: if the entry is not in `SVTYPE_contig:start-end` form
    """

    _svtype, _, locus = entry.partition('_')
    contig, _, span = locus.partition(':')
    start, _, end = span.partition('-')
    if not (contig and start and end):
        raise ValueError(f'Unparseable CPX_INTERVALS entry: {entry}')
    return contig, int(start), int(end)


def sorted_intervals(intervals: str) -> str:
    """
    Reorder one record's CPX_INTERVALS value into coordinate order.

    Args:
        intervals (str): the raw comma-separated CPX_INTERVALS value

    Returns:
        str: the same entries, coordinate-sorted, comma-separated

    Raises:
        ValueError: if any entry is not in `SVTYPE_contig:start-end` form
    """

    entries = intervals.split(',')

    contig_order: dict[str, int] = {}
    for entry in entries:
        contig, _start, _end = parse_interval(entry)
        contig_order.setdefault(contig, len(contig_order))

    return ','.join(sorted(entries, key=lambda entry: interval_sort_key(entry, contig_order)))


def cli_main():
    """
    main method wrapper for console script execution
    """
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, help='SV VCF to sort CPX_INTERVALS in')
    parser.add_argument('--output', required=True, help='Where to write the sorted VCF')
    args = parser.parse_args()
    main(vcf_in=args.input, vcf_out=args.output)


def main(vcf_in: str, vcf_out: str):
    """
    Coordinate-sort CPX_INTERVALS on every record that carries it, and write the VCF back out.

    Args:
        vcf_in (str): path to the SV VCF
        vcf_out (str): path to write the sorted VCF to
    """

    vcf = VCF(vcf_in)
    writer = Writer(vcf_out, vcf)

    variants = 0
    complex_variants = 0
    reordered = 0
    unparseable = 0

    for variant in vcf:
        variants += 1

        # cyvcf2 hands back the raw comma-joined string for this Number=. String field
        if intervals := variant.INFO.get(CPX_INTERVALS):
            complex_variants += 1
            try:
                if (in_order := sorted_intervals(intervals)) != intervals:
                    variant.INFO[CPX_INTERVALS] = in_order
                    reordered += 1
            except ValueError as e:
                # pass the record through untouched - SVAnnotate's own parser may still cope, and failing
                # the whole cohort over one unexpected entry format would be a worse outcome
                unparseable += 1
                logger.warning(f'{variant.ID}: leaving {CPX_INTERVALS} unsorted - {e}')

        writer.write_record(variant)

    writer.close()
    vcf.close()

    logger.info(f'Checked {complex_variants} complex variants across {variants} records')
    logger.info(f'Coordinate-sorted {CPX_INTERVALS} on {reordered} of them')

    if unparseable:
        logger.warning(f'{unparseable} records had an unparseable {CPX_INTERVALS} and were left untouched')


if __name__ == '__main__':
    cli_main()
