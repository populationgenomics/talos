"""
Hail Query functions for seqr loader.
"""

import logging
import os

import hail as hl

from cpg_utils.hail_batch import reference_path, genome_build


def annotate_cohort(
    vcf_path, out_mt_path, vep_ht_path, checkpoint_prefix=None, vep_only: bool = False
):
    """
    Convert VCF to matrix table, annotate for Seqr Loader, add VEP and VQSR annotations.

    Args:
        vcf_path ():
        out_mt_path ():
        vep_ht_path ():
        checkpoint_prefix ():
        vep_only ():

    Returns:

    """

    # set up a logger in this Hail Query runtime
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    def _read(path):
        if path.strip('/').endswith('.ht'):
            t = hl.read_table(str(path))
        else:
            assert path.strip('/').endswith('.mt')
            t = hl.read_matrix_table(str(path))
        logging.info(f'Read checkpoint {path}')
        return t

    def _checkpoint(t, file_name):
        if checkpoint_prefix:
            path = os.path.join(checkpoint_prefix, file_name)
            logging.info(f'Checkpointing to {path}')

            t.write(str(path), overwrite=True)
            logging.info(f'Wrote checkpoint {path}')
            t = _read(str(path))
        return t

    logger.info(f'Importing VCF {vcf_path}')
    mt = hl.import_vcf(
        str(vcf_path),
        reference_genome=genome_build(),
        skip_invalid_loci=True,
        force_bgz=True,
    )

    logger.info(f'Loading VEP Table from {vep_ht_path}')
    vep_ht = _read(vep_ht_path)

    logger.info(f'Adding VEP annotations into the Matrix Table from {vep_ht_path}')
    mt = mt.annotate_rows(vep=vep_ht[mt.locus, mt.alleles].vep)

    if vep_only:
        logger.info(f'Writing VEP-only annotation to {out_mt_path}')
        mt.describe()
        mt.write(str(out_mt_path), overwrite=True)
        return

    mt = _checkpoint(mt, 'mt-plus-vep.mt')

    # Add potentially missing fields
    if not all(attr in mt.row_value for attr in ['AC', 'AF', 'AN']):
        if mt.count_cols() == 0:
            logger.info('No samples in the Matrix Table, adding dummy values')
            mt = mt.annotate_rows(AN=1, AF=0.01, AC=1)
        else:
            logger.info('Adding AC/AF/AN attributes from variant_qc')
            mt = hl.variant_qc(mt)
            mt = mt.annotate_rows(
                AN=mt.variant_qc.AN, AF=mt.variant_qc.AF[1], AC=mt.variant_qc.AC[1]
            )
            mt = mt.drop('variant_qc')

    # don't fail if the AC/AF attributes are an inappropriate type
    for attr in ['AC', 'AF']:
        if isinstance(mt[attr], hl.ArrayExpression):
            mt = mt.annotate_rows(**{attr: mt[attr][1]})

    logger.info('Annotating with seqr-loader aggregate data')
    ref_ht = _read(str(reference_path('seqr_combined_reference_data')))
    clinvar_ht = _read(str(reference_path('seqr_clinvar')))
    mt = mt.annotate_rows(clinvar_data=clinvar_ht[mt.row_key], **ref_ht[mt.row_key])

    mt = _checkpoint(mt, 'mt-plus-ref-data.mt')

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

    logger.info(f'Written final matrix table into {out_mt_path}')
