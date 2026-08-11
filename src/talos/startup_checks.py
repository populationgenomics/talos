"""
This is a startup checker module

Takes the inputs (mandatory and optional) and checks if they are valid.

The small-variant input is a BCSQ + echtvar annotated VCF shard, and ClinVar arrives as the
ClinvArbitration release VCF, so the data checks here are header and first-rows probes read with
cyvcf2 - "taste the inputs" - rather than the MatrixTable schema walk this replaces.
"""

import sys
from argparse import ArgumentParser
from os import getenv

import pendulum
from cyvcf2 import VCF
from loguru import logger
from mendelbrot.pedigree_parser import PedigreeParser

from talos.config import config_check, config_retrieve
from talos.run_streaming_filtering import GNOMAD_SOURCE_FIELDS
from talos.vcf_streaming import header_has_field, split_csq_header

# collect all parsing errors as strings, print before crashing (unless everything passes...)
LOG_ERRORS: list[str] = []
CONFIG_ERRORS: list[str] = []

# INFO fields every annotated shard must declare, mapped to the process which should have written them
REQUIRED_INFO_FIELDS: dict[str, str] = {
    'AC': 'bcftools +fill-tags, in NormaliseVcf',
    'AN': 'bcftools +fill-tags, in NormaliseVcf',
    'BCSQ': 'bcftools csq, in AnnotateCsqWithBcftools',
    **dict.fromkeys(GNOMAD_SOURCE_FIELDS.values(), 'the gnomAD echtvar zip, in AnnotateWithEchtvar'),
}

# BCSQ fields read from every consequence - the rest of the csq string is derived after the split
REQUIRED_CSQ_FIELDS = ('consequence', 'gene', 'transcript', 'biotype', 'amino_acid_change', 'dna_change')

# fields structure_consequences adds post-split, so a configured csq_string can legitimately name them
DERIVED_CSQ_FIELDS = ('codon', 'am_class', 'am_pathogenicity', 'mane_status', 'ensp', 'mane_id', 'gene_id')

# INFO fields the ClinvArbitration release VCF must carry, aliased to clinvar_* by the echtvar encode
REQUIRED_CLINVAR_FIELDS = ('allele_id', 'gold_stars', 'clinical_significance')

# how many rows of the annotated VCF to read when probing BCSQ arity
CSQ_PROBE_ROWS = 500

# a ClinVar release with fewer rows than this is test data, not a release
CLINVAR_MIN_ROWS = 100


# Define expected types for each config key (update as needed)
SCHEMA = {
    'GeneratePanelData': {
        'default_panel': int,
        'panelapp': str,
    },
    'RunHailFiltering': {
        'ac_threshold': float,
        'additional_csq': list,
        'af_semi_rare': float,
        'callset_af_sv_recessive': float,
        'critical_csq': list,
        'minimum_depth': int,
        'csq_string': list,
        'de_novo': {
            'min_child_ab': float,
            'min_depth': int,
            'max_depth': int,
            'min_proband_gq': int,
            'min_alt_depth': int,
        },
    },
    'ValidateMOI': {
        'min_callset_ac_to_filter': int,
        'gnomad_max_af': float,
        'gnomad_sv_max_af': float,
        'callset_max_af': float,
        'callset_sv_max_af': float,
        'gnomad_max_homozygotes': int,
        'gnomad_max_hemizygotes': int,
        'dominant_gnomad_max_af': float,
        'dominant_gnomad_sv_max_af': float,
        'dominant_gnomad_max_ac': int,
        'dominant_gnomad_max_homozygotes': int,
        'dominant_callset_max_af': float,
        'dominant_callset_sv_max_af': float,
        'dominant_callset_max_ac': int,
        'clinvar_gnomad_max_af': float,
        'clinvar_dominant_gnomad_max_af': float,
        'clinvar_callset_max_af': float,
        'clinvar_dominant_callset_max_af': float,
        'ignore_categories': list,
        'support_categories': list,
        'phenotype_match': list,
    },
}
SCHEMA_OPTIONAL = {
    'GeneratePanelData': {
        'require_pheno_match': list,
        'forbidden_genes': list,
        'forced_panels': list,
        'manual_overrides': list,
        'within_x_months': int,
    },
    'RunHailFilteringSv': {
        'gnomad_population': str,
    },
    'ValidateMOI': {
        'ignore_categories': list,
        'support_categories': list,
        'phenotype_match': list,
    },
    'HPOFlagging': {
        'semantic_match': bool,
        'min_similarity': float,
    },
}


def check_csq_string_config(csq_fields: list[str]):
    """
    The output csq string is assembled from config, so a configured name which is neither a BCSQ
    field nor derived downstream renders as an empty string for every variant.

    This is a warning, not an error - an over-specified csq_string is cosmetic, not a data problem.

    Args:
        csq_fields (list[str]): field names parsed from this VCF's BCSQ header line
    """
    configured = config_retrieve(['RunHailFiltering', 'csq_string'], [])
    available = set(csq_fields) | set(DERIVED_CSQ_FIELDS)
    if unpopulated := [field for field in configured if field not in available]:
        logger.warning(f'csq_string names fields which are never populated, they will render empty: {unpopulated}')


def check_pedigree_overlap(reader: VCF, pedigree: PedigreeParser, vcf_path: str):
    """
    A run with no shared samples silently produces an empty report, so fail on it here instead.

    Args:
        reader (VCF): the opened annotated VCF
        pedigree (PedigreeParser): the parsed pedigree
        vcf_path (str): path to the VCF, for error messages
    """
    vcf_samples = set(reader.samples)
    if not vcf_samples & set(pedigree.get_all_sample_ids()):
        LOG_ERRORS.append(f'No pedigree samples are present in {vcf_path}')
    elif not vcf_samples & set(pedigree.get_affected_member_ids()):
        LOG_ERRORS.append(f'No affected pedigree samples are present in {vcf_path}')


def probe_csq_arity(reader: VCF, csq_fields: list[str], vcf_path: str):
    """
    Read the first rows of the VCF and check the BCSQ payload against its own header.

    BCFtools drops trailing empty fields, so a short entry is expected, but an entry with more
    fields than the header describes means data and header have drifted apart, and every
    consequence would be parsed against the wrong field names.

    Args:
        reader (VCF): the opened annotated VCF, positioned at the first row
        csq_fields (list[str]): field names parsed from this VCF's BCSQ header line
        vcf_path (str): path to the VCF, for error messages
    """
    rows = 0
    rows_with_bcsq = 0

    for variant in reader:
        rows += 1
        if bcsq := variant.INFO.get('BCSQ'):
            rows_with_bcsq += 1
            for entry in bcsq.split(','):
                if len(entry.split('|')) > len(csq_fields):
                    LOG_ERRORS.append(
                        f'BCSQ at {variant.CHROM}:{variant.POS} has more fields than the '
                        f'{len(csq_fields)} its header declares: {entry}',
                    )
                    break
        if rows >= CSQ_PROBE_ROWS:
            break

    if not rows:
        LOG_ERRORS.append(f'VCF contains no variants: {vcf_path}')
    elif not rows_with_bcsq:
        LOG_ERRORS.append(f'None of the first {rows} variants in {vcf_path} carry a BCSQ annotation')


def check_vcf(vcf_path: str, pedigree: PedigreeParser | None):
    """
    Check the annotated VCF declares every INFO field the streaming filter reads, and that the
    BCSQ contents match the field list its header advertises.

    One shard is representative - every shard is written by the same processes, over the same
    samples, so a contract failure here is a contract failure everywhere.

    Args:
        vcf_path (str): path to one annotated VCF shard
        pedigree (PedigreeParser | None): the parsed pedigree, if it was readable
    """
    logger.info(f'Checking annotated VCF {vcf_path}')

    try:
        reader = VCF(vcf_path)
    except OSError as ose:
        LOG_ERRORS.append(f'VCF could not be opened: {vcf_path}\n{ose}')
        return

    for field, source in REQUIRED_INFO_FIELDS.items():
        if not header_has_field(reader, field):
            LOG_ERRORS.append(f'INFO/{field} is missing from {vcf_path}, expected from {source}')

    if pedigree is not None:
        check_pedigree_overlap(reader, pedigree, vcf_path)

    # without a BCSQ header line there is no field list to check the rows against
    if not header_has_field(reader, 'BCSQ'):
        return

    csq_fields = split_csq_header(reader)
    if missing := [field for field in REQUIRED_CSQ_FIELDS if field not in csq_fields]:
        LOG_ERRORS.append(f'BCSQ header of {vcf_path} does not describe required fields: {missing}')

    check_csq_string_config(csq_fields)
    probe_csq_arity(reader, csq_fields, vcf_path)


def clinvar_file_date(reader: VCF) -> str | None:
    """Pull the ##fileDate value out of a raw VCF header, if it carries one."""
    for line in reader.raw_header.splitlines():
        if line.startswith('##fileDate='):
            return line.removeprefix('##fileDate=').strip()
    return None


def check_clinvar_freshness(reader: VCF, clinvar_path: str):
    """
    Check the ClinVar release is recent, reading ##fileDate from the header. This replaces the
    creation_date global of the Hail Table the release VCF supersedes.

    Args:
        reader (VCF): the opened ClinvArbitration release VCF
        clinvar_path (str): path to the VCF, for error messages
    """
    if (file_date := clinvar_file_date(reader)) is None:
        LOG_ERRORS.append(f'ClinVar VCF lacks a ##fileDate header line: {clinvar_path}')
        return

    try:
        created = pendulum.from_format(file_date, 'YYYY-MM-DD')
    except ValueError:
        LOG_ERRORS.append(f'ClinVar VCF ##fileDate is not a YYYY-MM-DD date: {file_date}')
        return

    if created < pendulum.now().subtract(months=2) and config_retrieve(['clinvar_check_age'], True):
        LOG_ERRORS.append(
            f'ClinVar VCF {clinvar_path} is > 2 months old: {created.to_date_string()}, get a new one.'
            f'Alternatively, disable this check by setting the config key "clinvar_check_age" to False.',
        )


def check_clinvar(clinvar_path: str | None):
    """
    Check the ClinvArbitration release VCF is readable, carries the INFO fields the echtvar encode
    reads, holds more than a test-data volume of decisions, and is recent.

    Args:
        clinvar_path (str | None): path to the ClinvArbitration release VCF
    """
    if clinvar_path is None:
        LOG_ERRORS.append('ClinVar path is not provided.')
        return

    logger.info(f'Checking ClinVar VCF {clinvar_path}')

    try:
        reader = VCF(clinvar_path)
    except OSError as ose:
        LOG_ERRORS.append(f'ClinVar VCF could not be opened: {clinvar_path}\n{ose}')
        return

    for field in REQUIRED_CLINVAR_FIELDS:
        if not header_has_field(reader, field):
            LOG_ERRORS.append(f'ClinVar VCF {clinvar_path} does not declare the INFO field {field}')

    check_clinvar_freshness(reader, clinvar_path)

    # stop counting once the file has proven itself - a real release is millions of rows
    rows = 0
    for _variant in reader:
        rows += 1
        if rows >= CLINVAR_MIN_ROWS:
            break
    if rows < CLINVAR_MIN_ROWS:
        LOG_ERRORS.append(f'ClinVar VCF has fewer than {CLINVAR_MIN_ROWS} entries: {clinvar_path}')


def validate_types(config, schema, path=''):
    errors = []
    for key, expected_type in schema.items():
        if key not in config:
            continue
        value = config[key]
        if isinstance(expected_type, dict):
            errors += validate_types(value, expected_type, path + key + '.')
        elif isinstance(expected_type, list):
            if not isinstance(value, list):
                errors.append(f'{path}{key}: expected list, got {type(value).__name__}')
        elif not isinstance(value, expected_type):
            errors.append(f'{path}{key}: expected {expected_type.__name__}, got {type(value).__name__}')
    return errors


def validate_pedigree(pedigree_path: str | None) -> PedigreeParser | None:
    """
    Check that the file exists, is readable, and contains affected members.

    Returns:
        the parsed pedigree, or None if it could not be read
    """
    if pedigree_path is None:
        LOG_ERRORS.append('Pedigree path is not provided.')
        return None

    logger.info(f'Checking pedigree at {pedigree_path}')

    # check if the pedigree file is readable
    try:
        pedigree = PedigreeParser(pedigree_path)
    except (ValueError, OSError) as error:
        LOG_ERRORS.append(f'Error parsing pedigree file: {pedigree_path}\n{error}')
        return None

    if not pedigree.get_affected_member_ids():
        LOG_ERRORS.append(f'Pedigree file is empty or does not contain affected members: {pedigree_path}')

    return pedigree


def recursive_schema_validation(schema: dict, lead: list[str] | None = None, optional: bool = False):
    """
    Validate the schema against the provided configuration. Re-call this method at each level of the schema.

    Args:
        schema (dict): the schema to use for validation
        lead (list[str] | None): the path to the current schema level, used for error messages and config retrieval
        optional (bool): if True, we will not raise an error if the key is missing, but will check the types if present
    """
    # this is used to build the keys
    if lead is None:
        lead = []

    for key, value in schema.items():
        if key not in config_retrieve(lead):
            if not optional:
                LOG_ERRORS.append(f'Missing required config key: {".".join([*lead, key])}')
            continue
        if isinstance(value, dict):
            # if the value is a dict, we need to recurse into it
            recursive_schema_validation(schema=value, lead=[*lead, key], optional=optional)
            continue

        key_path = [*lead, key]

        CONFIG_ERRORS.extend(config_check(key=key_path, expected_type=value, optional=optional))


def check_config():
    """
    Configuration checks, we fire off a method to run recursively through the schema and check the types.
    If the environment variable TALOS_CONFIG is not set, we will not be able to access, so it's a hard fail.
    """

    if (config_path := getenv('TALOS_CONFIG')) is None:
        LOG_ERRORS.append("Environment variable TALOS_CONFIG is not set, config won't be accessible.")
        return

    logger.info(f'Checking config at {config_path}')

    recursive_schema_validation(SCHEMA, optional=False)
    recursive_schema_validation(SCHEMA_OPTIONAL, optional=True)

    # add the config-specific errors to the logging errors
    LOG_ERRORS.extend(CONFIG_ERRORS)


def main(
    pedigree_path: str | None,
    vcf_path: str,
    clinvar_path: str | None,
) -> None:
    """
    Main function to run all startup checks.
    """

    # Run the checks on the pedigree, and on the Config file (picked up from the environment variable TALOS_CONFIG)
    pedigree = validate_pedigree(pedigree_path)
    check_config()

    check_vcf(vcf_path, pedigree)
    check_clinvar(clinvar_path)

    if LOG_ERRORS:
        logger.error('One or more startup checks failed:')
        logger.error('\n'.join(LOG_ERRORS))
        sys.exit(1)

    print('All startup checks passed successfully.')


if __name__ == '__main__':
    parser = ArgumentParser(description='Startup checks for Talos pipeline')
    parser.add_argument('--pedigree', help='Path to the pedigree file.', default=None)
    parser.add_argument('--vcf', help='Path to one annotated VCF shard.', required=True)
    parser.add_argument('--clinvar', help='Path to the ClinvArbitration release VCF.', default=None)
    args = parser.parse_args()
    main(pedigree_path=args.pedigree, vcf_path=args.vcf, clinvar_path=args.clinvar)
