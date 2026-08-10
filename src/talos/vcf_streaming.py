"""
Shared helpers for the cyvcf2 streaming pipeline stages.

These replace the Hail-based reshaping in annotated_vcf_into_matrixtable.py and
run_hail_filtering.py with plain-Python equivalents operating on one variant at
a time. See docs/proposals/reducing-spark-dependence.md.
"""

import gzip
import re
from typing import TYPE_CHECKING, Any

from loguru import logger
from mendelbrot.pedigree_parser import PedigreeParser

from talos.config import config_retrieve

if TYPE_CHECKING:
    import cyvcf2

# ClinvArbitration decision strings (see run_hail_filtering.py constants)
BENIGN = 'benign'
PATHOGENIC = 'Pathogenic/Likely Pathogenic'

# sentinel used across the labelled VCFs for absent String values
MISSING_STRING = 'missing'

# BCSQ fields are pipe-delimited; a fully populated entry has 7 fields, but
# BCFtools truncates the string when trailing fields are absent
BCSQ_EXPECTED_FIELDS = 7

LEADING_INT_RE = re.compile(r'^(\d+)')

CLINVAR_TSV_KEYS = ['contig', 'position', 'reference', 'alternate', 'clinical_significance', 'gold_stars', 'allele_id']


def normalise_chrom(chrom: str) -> str:
    """Reduce a contig name to a chr-less, MT-as-M form for cross-source key matching."""
    chrom = chrom.removeprefix('chr')
    return 'M' if chrom == 'MT' else chrom


def split_csq_header(reader: 'cyvcf2.VCF') -> list[str]:
    """
    Pull the BCSQ header line from an opened VCF and split the Format description
    into a list of lowercase field names.

    Equivalent to annotated_vcf_into_matrixtable.extract_and_split_csq_string.
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
        csq['codon'] = int(codon_match.group(1)) if codon_match else None

        consequences.append(csq)

    return consequences


def consequences_to_csq_string(consequences: list[dict[str, Any]]) -> str:
    """
    Collapse consequence dicts into the single INFO-ready csq String.

    Field selection and ordering come from config (RunHailFiltering.csq_string),
    matching run_hail_filtering.csq_struct_to_string. Absent/None values render
    as empty strings.
    """

    csq_fields = config_retrieve(['RunHailFiltering', 'csq_string'])

    entries = []
    for csq in consequences:
        values = [csq.get(field) for field in csq_fields]
        entries.append('|'.join('' if value is None else str(value) for value in values))
    return ','.join(entries)


def read_clinvar_decisions_tsv(tsv_path: str, contigs: set[str] | None = None) -> dict[tuple, dict[str, Any]]:
    """
    Read a ClinvArbitration decisions TSV into a lookup dict.

    Columns: contig, position, reference, alternate, clinical_significance, gold_stars, allele_id.

    Args:
        tsv_path (str): path to the (optionally gzipped) decisions TSV
        contigs (set[str] | None): if provided, only keep entries on these
            normalised contigs (see normalise_chrom)

    Returns:
        {(contig, pos, ref, alt): {'clinical_significance': str, 'gold_stars': int, 'allele_id': int}}
    """

    decisions: dict[tuple, dict[str, Any]] = {}

    open_method = gzip.open if tsv_path.endswith('.gz') else open
    with open_method(tsv_path, 'rt') as handle:
        header = handle.readline().rstrip('\n').split('\t')
        if unknown := set(header) - set(CLINVAR_TSV_KEYS):
            raise ValueError(f'Unexpected columns in ClinVar decisions TSV {tsv_path}: {unknown}')

        for line in handle:
            entry = dict(zip(header, line.rstrip('\n').split('\t'), strict=True))
            contig = normalise_chrom(entry['contig'])
            if contigs and contig not in contigs:
                continue
            key = (contig, int(entry['position']), entry['reference'], entry['alternate'])
            decisions[key] = {
                'clinical_significance': entry['clinical_significance'],
                'gold_stars': int(entry['gold_stars']),
                'allele_id': int(entry['allele_id']),
            }

    logger.info(f'Read {len(decisions)} ClinVar decisions from {tsv_path}')
    return decisions


def variant_is_pass(variant: 'cyvcf2.Variant') -> bool:
    """cyvcf2 represents PASS/empty FILTER as None."""
    return variant.FILTER is None


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
