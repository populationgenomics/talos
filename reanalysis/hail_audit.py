"""
pulled out the methods relating to the 'hail audit'
mostly just to reduce the line count
partially to make this available in other places
"""


import logging

import hail as hl

BASE_FIELDS_REQUIRED = [
    ('locus', hl.LocusExpression),
    ('alleles', hl.ArrayExpression),
    ('AC', hl.Int32Expression),
    ('AF', hl.Float64Expression),
    ('AN', hl.Int32Expression),
]
FIELDS_REQUIRED = {
    'splice_ai': [
        ('delta_score', hl.Float32Expression),
        ('splice_consequence', hl.StringExpression),
    ],
    'gnomad_exomes': [
        ('AF', hl.Float64Expression),
        ('AN', hl.Int32Expression),
        ('AC', hl.Int32Expression),
        ('Hom', hl.Int32Expression),
        ('Hemi', hl.Int32Expression),
    ],
    'gnomad_genomes': [
        ('AF', hl.Float64Expression),
        ('AN', hl.Int32Expression),
        ('AC', hl.Int32Expression),
        ('Hom', hl.Int32Expression),
        ('Hemi', hl.Int32Expression),
    ],
    'cadd': [('PHRED', hl.Float32Expression)],
    'dbnsfp': [
        ('REVEL_score', hl.StringExpression),
        ('MutationTaster_pred', hl.StringExpression),
    ],
    'clinvar': [
        ('clinical_significance', hl.StringExpression),
        ('gold_stars', hl.Int32Expression),
    ],
    'geneIds': [],
}

VEP_TX_FIELDS_REQUIRED = [
    ('variant_allele', hl.StringExpression),
    ('consequence_terms', hl.ArrayExpression),
    ('transcript_id', hl.StringExpression),
    ('protein_id', hl.StringExpression),
    ('gene_id', hl.StringExpression),
    ('gene_symbol', hl.StringExpression),
    ('gene_symbol_source', hl.StringExpression),
    ('canonical', hl.Int32Expression),
    ('cdna_start', hl.Int32Expression),
    ('cds_start', hl.Int32Expression),
    ('cds_end', hl.Int32Expression),
    ('biotype', hl.StringExpression),
    ('protein_start', hl.Int32Expression),
    ('protein_end', hl.Int32Expression),
    ('sift_score', hl.Float64Expression),
    ('sift_prediction', hl.StringExpression),
    ('polyphen_score', hl.Float64Expression),
    ('mane_select', hl.StringExpression),
    ('lof', hl.StringExpression),
]

USELESS_FIELDS = [
    'a_index',
    'aIndex',
    'alleles_old',
    'AS_culprit',
    'clinvar_data',
    'docId',
    'domains',
    'eigen',
    'g1k',
    'geno2mp',
    'genotypes',
    'gnomad_exome_coverage',
    'gnomad_genome_coverage',
    'locus_old',
    'mainTranscript',
    'mpc',
    'negative_train_site',
    'originalAltAlleles',
    'positive_train_site',
    'primate_ai',
    'rg37_locus',
    'samples_ab',
    'samples_gq',
    'samples_no_call',
    'samples_num_alt',
    'score',
    'sortedTranscriptConsequences',
    'topmed',
    'transcriptConsequenceTerms',
    'transcriptIds',
    'was_split',
    'wasSplit',
    'vep_proc_id',
]


def fields_audit(
    mt: hl.MatrixTable, base_fields: list[tuple], nested_fields: dict
) -> bool:
    """
    checks that the required fields are all present before continuing
    """
    problems = []
    # iterate over top-level attributes
    for annotation, datatype in base_fields:
        if annotation in mt.row_value or annotation in mt.row_key:
            if not isinstance(mt[annotation], datatype):
                problems.append(f'{annotation}: {datatype}/{type(mt[annotation])}')
        else:
            problems.append(f'{annotation}:missing')

    for field_group, group_types in nested_fields.items():
        if field_group not in mt.row_value:
            problems.append(f'{field_group}:missing')
        else:
            for annotation, datatype in group_types:
                if annotation in mt[field_group]:
                    if not isinstance(mt[field_group][annotation], datatype):
                        problems.append(
                            f'{annotation}:'
                            f'{datatype}/'
                            f'{type(mt[field_group][annotation])}'
                        )
                else:
                    problems.append(f'{annotation}:missing')
    for problem in problems:
        logging.error(f'MT field: \t{problem}')
    return len(problems) == 0


def vep_audit(
    mt: hl.MatrixTable, expected_fields: list[tuple[str, hl.Expression]]
) -> bool:
    """
    check that the required VEP annotations are present
    True if the 'audit' passes (all required fields present)
    """

    problems = []
    # now the content of the transcript_consequences
    if 'vep' not in mt.row_value:
        problems.append('VEP:missing')
    elif 'transcript_consequences' not in mt.vep:
        problems.append('transcript_consequences:missing')
    else:
        fields_and_types = dict(mt.vep.transcript_consequences[0].items())
        for field, field_type in expected_fields:
            if field in fields_and_types:
                if not isinstance(fields_and_types[field], field_type):
                    problems.append(
                        f'{field}:{field_type}/{type(fields_and_types[field])}'
                    )
            else:
                problems.append(f'{field}:missing')

    for problem in problems:
        logging.error(f'VEP field: {problem}')

    return len(problems) == 0
