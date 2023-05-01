"""
Hail Query functions for seqr loader.
"""

import logging
import os

import hail as hl

from cpg_utils.hail_batch import reference_path, genome_build


def annotate_cohort(
    vcf_path,
    out_mt_path,
    vep_ht_path,
    checkpoint_prefix=None,
):
    """
    Convert VCF to matrix table, annotate for Seqr Loader, add VEP and VQSR
    annotations.
    """

    mt = hl.import_vcf(
        str(vcf_path),
        reference_genome=genome_build(),
        skip_invalid_loci=True,
        force_bgz=True,
    )
    logging.info(f'Importing VCF {vcf_path}')
    if checkpoint_prefix:
        mt = mt.checkpoint(
            os.path.join(checkpoint_prefix, 'mt-from-vcf.mt'), overwrite=True
        )

    logging.info(f'Loading VEP Table from {vep_ht_path}')
    vep_ht = hl.read_table(vep_ht_path)
    logging.info(f'Adding VEP annotations into the Matrix Table from {vep_ht_path}')
    mt = mt.annotate_rows(vep=vep_ht[mt.locus, mt.alleles].vep)

    # Add potentially missing fields
    if not all(attr in mt.row_value for attr in ['AC', 'AF', 'AN']):
        if mt.count_cols() == 0:
            logging.info('No samples in the Matrix Table, adding dummy values')
            mt = mt.annotate_rows(AN=1, AF=0.01, AC=1)
        else:
            logging.info('Adding AC/AF/AN attributes from variant_qc')
            mt = hl.variant_qc(mt)
            mt = mt.annotate_rows(
                AN=mt.variant_qc.AN, AF=mt.variant_qc.AF[1], AC=mt.variant_qc.AC[1]
            )
            mt = mt.drop('variant_qc')

    # don't fail if the AC/AF attributes are an inappropriate type
    for attr in ['AC', 'AF']:
        if isinstance(mt[attr], hl.ArrayExpression):
            mt = mt.annotate_rows(**{attr: mt[attr][1]})

    logging.info('Annotating with seqr-loader aggregate data')
    ref_ht = hl.read_table(str(reference_path('seqr_combined_reference_data')))
    clinvar_ht = hl.read_table(str(reference_path('seqr_clinvar')))
    mt = mt.annotate_rows(clinvar_data=clinvar_ht[mt.row_key], **ref_ht[mt.row_key])

    mt = mt.annotate_rows(
        geneIds=hl.set(mt.vep.transcript_consequences.map(lambda c: c.gene_id)),
        clinvar=hl.struct(
            **{
                'allele_id': mt.clinvar_data.info.ALLELEID,
                'clinical_significance': hl.delimit(mt.clinvar_data.info.CLNSIG),
                'gold_stars': mt.clinvar_data.gold_stars,
            }
        ),
    )

    mt = mt.annotate_globals(
        sourceFilePath=vcf_path,
        genomeVersion=genome_build().replace('GRCh', ''),
        hail_version=hl.version(),
    )

    mt.describe()
    mt.write(str(out_mt_path), overwrite=True)

    logging.info(f'Written final matrix table into {out_mt_path}')
