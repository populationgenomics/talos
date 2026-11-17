"""
Shared helpers for the cyvcf2 streaming pipeline stages.
"""

import re
from typing import TYPE_CHECKING, Any

from loguru import logger
from mendelbrot.pedigree_parser import PedigreeParser

from talos.config import config_retrieve

if TYPE_CHECKING:
    import cyvcf2

# ClinvArbitration decision strings.
# ClinvArbitration 3.0.0 shortened its P/LP rating from 'Pathogenic/Likely Pathogenic' to
# 'Pathogenic', so is_pathogenic matches the shared word rather than comparing for equality -
# an equality check silently stops labelling anything the next time the wording moves
BENIGN = 'benign'
PATHOGENIC = 'Pathogenic'

# sentinel used across the labelled VCFs for absent String values
MISSING_STRING = 'missing'

# BCSQ fields are pipe-delimited; a fully populated entry has 7 fields, but
# BCFtools truncates the string when trailing fields are absent
BCSQ_EXPECTED_FIELDS = 7

LEADING_INT_RE = re.compile(r'^(\d+)')


def normalise_chrom(chrom: str) -> str:
    """Reduce a contig name to a chr-less, MT-as-M form for cross-source key matching."""
    chrom = chrom.removeprefix('chr')
    return 'M' if chrom == 'MT' else chrom


def split_csq_header(reader: 'cyvcf2.VCF') -> list[str]:
    """
    Pull the BCSQ header line from an opened VCF and split the Format description into a list of lowercase field names.
    """
    bcsq_header = reader.get_header_type('BCSQ')
    description = bcsq_header['Description'].strip('"')
    return description.split('Format: ')[-1].lower().split('|')


def parse_bcsq_entries(bcsq_string: str, csq_fields: list[str]) -> list[dict[str, Any]]:
    """
    Split a raw BCSQ INFO value into one dict per consequence.

    - pads truncated entries (BCFtools drops trailing empty fields)
    - drops the strand field
    - derives an integer codon from amino_acid_change ("123P" or "123P-124F" -> 123)

    Args:
        bcsq_string (str): the raw comma-joined BCSQ INFO value
        csq_fields (list[str]): lowercase field names from the BCSQ header

    Returns:
        list of consequence dicts, keyed by BCSQ field name plus 'codon'
    """

    consequences = []
    for entry in bcsq_string.split(','):
        values = entry.split('|')
        # pad truncated entries out to the full field count
        values.extend([''] * (len(csq_fields) - len(values)))

        csq = {field: value for field, value in zip(csq_fields, values, strict=False) if field != 'strand'}

        codon_match = LEADING_INT_RE.match(csq.get('amino_acid_change', ''))
        csq['codon'] = int(codon_match.group(1)) if codon_match else None  # type: ignore[assignment]

        consequences.append(csq)

    return consequences


def consequences_to_csq_string(consequences: list[dict[str, Any]]) -> str:
    """
    Collapse consequence dicts into the single INFO-ready csq String.

    Field selection and ordering come from config (RunSmallFiltering.csq_string).
    Absent/None values render as empty strings.
    """

    csq_fields = config_retrieve(['RunSmallFiltering', 'csq_string'])

    entries = []
    for csq in consequences:
        values = [csq.get(field) for field in csq_fields]
        entries.append('|'.join('' if value is None else str(value) for value in values))
    return ','.join(entries)


def read_clinvar_fields(info: Any) -> tuple[str, int, int]:
    """
    Normalise the echtvar-applied ClinvArbitration INFO fields.
    Absent annotations (or echtvar missing-value sentinels) become 'missing'/0.
    """
    significance = info.get('clinvar_significance')
    if significance in (None, '', '.'):
        significance = MISSING_STRING
    # VCF INFO values cannot contain spaces, so the echtvar zip encodes them as underscores
    significance = significance.replace('_', ' ')
    stars = info.get('clinvar_stars')
    stars = int(stars) if stars is not None and stars >= 0 else 0
    allele = info.get('clinvar_allele')
    allele = int(allele) if allele is not None and allele >= 0 else 0
    return significance, stars, allele


def is_pathogenic(significance: str) -> bool:
    """
    Is this ClinvArbitration rating a Pathogenic/Likely-Pathogenic one?

    Matched as a substring so that both the 2.x rating ('Pathogenic/Likely Pathogenic') and the
    3.0 rating ('Pathogenic') are recognised. No other rating - Benign, Conflicting, VUS,
    Unknown - contains the word.
    """
    return PATHOGENIC.lower() in significance.lower()


def variant_is_pass(variant: 'cyvcf2.Variant') -> bool:
    """cyvcf2 represents PASS/empty FILTER as None."""
    return variant.FILTER is None or (variant.INFO.get('SVTYPE') == 'BND' and variant.FILTERS == ['UNRESOLVED'])


def first_value(value: Any) -> Any:
    """cyvcf2 collapses single-element Number=A INFO fields to a scalar - undo that."""
    return value[0] if isinstance(value, tuple | list) else value


def header_has_field(reader: 'cyvcf2.VCF', field_id: str) -> bool:
    """Check whether the opened VCF declares a header entry with this ID."""
    try:
        reader.get_header_type(field_id)
    except KeyError:
        return False
    return True


def parse_pedigree(pedigree_path: str) -> PedigreeParser:
    """Parse the pedigree, reducing to affected singletons if the config says so."""
    pedigree_data = PedigreeParser(pedigree_path)
    if config_retrieve('singletons', False):
        logger.info('Reducing pedigree to affected singletons only')
        pedigree_data.set_participants(pedigree_data.as_singletons())
        pedigree_data.set_participants(pedigree_data.get_affected_members())
    return pedigree_data


def subset_reader_to_pedigree(reader: 'cyvcf2.VCF', pedigree_data: PedigreeParser) -> set[str]:
    """
    Restrict an opened VCF reader to the samples shared with the pedigree.

    Returns the set of shared samples - empty when there is no overlap, in which
    case the reader is left untouched and the caller decides how to bail out.
    """
    vcf_samples = set(reader.samples)
    ped_samples = pedigree_data.get_all_sample_ids()
    common_samples = ped_samples & vcf_samples
    logger.info(f'Samples in Pedigree: {len(ped_samples)}, VCF: {len(vcf_samples)}, common: {len(common_samples)}')
    if common_samples and common_samples != vcf_samples:
        reader.set_samples(sorted(common_samples))
    return common_samples
