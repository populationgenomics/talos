"""
This is a startup checker module

Takes the inputs (mandatory and optional) and checks if they are valid.
"""

import sys
from argparse import ArgumentParser
from os import getenv

import pendulum
from cloudpathlib.anypath import to_anypath
from loguru import logger
from mendelbrot.pedigree_parser import PedigreeParser

import hail as hl

from talos.config import config_check, config_retrieve

# collect all parsing errors as strings, print before crashing (unless everything passes...)
LOG_ERRORS: list[str] = []
CONFIG_ERRORS: list[str] = []

# hierarchy and types of the fields we expect in the MatrixTable
BASE_FIELDS_REQUIRED = [
    ('locus', hl.LocusExpression),
    ('alleles', hl.ArrayExpression),
    ('gene_ids', hl.SetExpression),
]
NESTED_FIELDS_REQUIRED = {
    'gnomad': [
        ('gnomad_AC', hl.Int32Expression),
        ('gnomad_AF', hl.Float64Expression),
        ('gnomad_HomAlt', hl.Int32Expression),
    ],
    'info': [
        ('AC', hl.ArrayNumericExpression),
        ('AF', hl.ArrayNumericExpression),
        ('AN', hl.Int32Expression),
    ],
    'transcript_consequences': [
        ('consequence', hl.StringExpression),
        ('gene', hl.StringExpression),
        ('gene_id', hl.StringExpression),
        ('transcript', hl.StringExpression),
        ('biotype', hl.StringExpression),
        ('amino_acid_change', hl.StringExpression),
        ('dna_change', hl.StringExpression),
        ('codon', hl.Int32Expression),
        ('am_class', hl.StringExpression),
        ('am_pathogenicity', hl.Float64Expression),
        ('mane_status', hl.StringExpression),
        ('ensp', hl.StringExpression),
        ('mane_id', hl.StringExpression),
    ],
}


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


def check_mt(mt_path: str | None):
    """
    Check that the MatrixTable exists and is readable.
    Then check the schema and data type of each element, this is a manual validation and kinda gross.
    """
    if mt_path is None:
        LOG_ERRORS.append('MatrixTable path is not provided.')
        return

    logger.info(f'Checking MatrixTable at {mt_path}')

    mt_anypath = to_anypath(mt_path)
    if not (mt_anypath / '_SUCCESS').exists():
        LOG_ERRORS.append(f'MatrixTable Success file does not exist: {mt_path}')
        return

    # boot it up and check the contents
    mt = hl.read_matrix_table(mt_path)

    # check the base fields - this works, but is only the key
    for field, field_type in BASE_FIELDS_REQUIRED:
        if field in mt.row_key or field in mt.row_value:
            if not isinstance(mt[field], field_type):
                LOG_ERRORS.append(f'{field}:{field_type}/{type(mt[field])}')
        else:
            LOG_ERRORS.append(f'{field}:missing from MT key')

    for field_group, group_types in NESTED_FIELDS_REQUIRED.items():
        if field_group not in mt.row_value:
            LOG_ERRORS.append(f'{field_group}:missing')
        else:
            # grab the first element from this array of structs as representative
            if field_group == 'transcript_consequences':
                fields_and_types = dict(mt.transcript_consequences[0].items())
            else:
                fields_and_types = dict(mt[field_group].items())
            for annotation, datatype in group_types:
                if annotation in fields_and_types:
                    if not isinstance(fields_and_types[annotation], datatype):
                        LOG_ERRORS.append(
                            f'{field_group}.{annotation}: {datatype}/{type(fields_and_types[annotation])}',
                        )
                else:
                    LOG_ERRORS.append(f'{annotation}:missing')


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


def validate_pedigree(pedigree_path: str | None):
    """
    Check that the file exists, is readable, and contains affected members.
    """
    if pedigree_path is None:
        LOG_ERRORS.append('Pedigree path is not provided.')
        return

    logger.info(f'Checking pedigree at {pedigree_path}')

    ped_path = to_anypath(pedigree_path)

    # quick existence check
    if not ped_path.exists():
        LOG_ERRORS.append(f'Pedigree file does not exist: {ped_path}')

    # check if the pedigree file is readable
    try:
        pedigree = PedigreeParser(pedigree_path)
    except ValueError as ve:
        LOG_ERRORS.append(f'Error parsing pedigree file: {ped_path}\n{ve}')
        return

    if not pedigree.get_affected_member_ids():
        LOG_ERRORS.append(f'Pedigree file is empty or does not contain affected members: {ped_path}')


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


def check_clinvar(clinvar_paths: list[str] | None):
    """
    Check that the ClinVar HailTable exists and is readable.
    """
    if clinvar_paths is None or (isinstance(clinvar_paths, list) and len(clinvar_paths) == 0):
        LOG_ERRORS.append('ClinVar paths are not provided.')
        return

    for clinvar_path in clinvar_paths:
        clinvar_anypath = to_anypath(clinvar_path)
        if not (clinvar_anypath / '_SUCCESS').exists():
            LOG_ERRORS.append(f'ClinVar Success file does not exist: {clinvar_path}')
            continue

        # boot it up and check the contents - we're looking for a non-empty table, much larger than the test-data
        clinvar_ht = hl.read_table(clinvar_path)
        if clinvar_ht.count() < 100:
            LOG_ERRORS.append(f'ClinVar HailTable has fewer than 100 entries: {clinvar_path}')
        if 'creation_date' not in clinvar_ht.globals:
            LOG_ERRORS.append('ClinVar HailTable lacks a creation date')
            continue
        created = pendulum.from_format(hl.eval(clinvar_ht.globals.creation_date), 'YYYY-MM-DD')
        if created < pendulum.now().subtract(months=2) and config_retrieve(['clinvar_check_age'], True):
            LOG_ERRORS.append(
                f'ClinVar HailTable {clinvar_path} is > 2 months old: {created.to_date_string()}, get a new one.'
                f'Alternatively, disable this check by setting the config key "clinvar_check_age" to False.',
            )


def main(
    pedigree_path: str | None,
    mt_path: str | None,
    clinvar_paths: list[str] | None,
) -> None:
    """
    Main function to run all startup checks.
    """

    # start up a hail runtime to check the MatrixTable and ClinVar HailTable
    hl.context.init_spark(master='local[*]', default_reference='GRCh38', quiet=True)

    # Run the checks on the pedigree, and on the Config file (picked up from the environment variable TALOS_CONFIG)
    validate_pedigree(pedigree_path)
    check_config()

    check_mt(mt_path)
    check_clinvar(clinvar_paths)

    if LOG_ERRORS:
        logger.error('One or more startup checks failed:')
        logger.error('\n'.join(LOG_ERRORS))
        sys.exit(1)

    print('All startup checks passed successfully.')


if __name__ == '__main__':
    parser = ArgumentParser(description='Startup checks for Talos pipeline')
    parser.add_argument('--pedigree', help='Path to the pedigree file.', default=None)
    parser.add_argument('--mt', help='Path to the MatrixTable.', default=None)
    parser.add_argument('--clinvar', help='Path to the ClinVar HailTable.', default=None, nargs='+')
    args = parser.parse_args()
    main(pedigree_path=args.pedigree, mt_path=args.mt, clinvar_paths=args.clinvar)
