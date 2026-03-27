"""
Merge pipeline ResultData with analyst decisions from the database.
"""

from talos.models import ReportVariant, ResultData


def variant_key(rv: ReportVariant) -> tuple:
    """Build the natural key for a variant, matching the decisions table unique constraint."""
    coords = rv.var_data.coordinates
    return (rv.sample, coords.chrom, coords.pos, coords.ref, coords.alt, rv.gene)


def collect_variant_keys(result_data: ResultData) -> list[tuple]:
    """Extract all variant natural keys from a ResultData object."""
    keys = []
    for participant in result_data.results.values():
        for rv in participant.variants:
            keys.append(variant_key(rv))
    return keys


def join_variants_with_decisions(
    result_data: ResultData,
    decisions: dict[tuple, dict],
) -> dict[str, list[dict]]:
    """
    For each sample in the result, attach matching decisions to their variants.

    Returns a dict of sample_id -> list of {report_variant, decision} dicts,
    preserving the structure needed for template rendering.
    """
    joined: dict[str, list[dict]] = {}

    for sample_id, participant in result_data.results.items():
        sample_variants = []
        for rv in participant.variants:
            key = variant_key(rv)
            sample_variants.append(
                {
                    'report_variant': rv,
                    'decision': decisions.get(key),
                }
            )
        joined[sample_id] = sample_variants

    return joined
