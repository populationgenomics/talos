"""
method file for re-sorting clinvar annotations by codon
utilises the latest in-house clinvar summary data
"""
import logging

import hail as hl

from reanalysis.hail_filter_and_label import PATHOGENIC


def protein_indexed_clinvar(mt: hl.MatrixTable, write_path: str, new_decisions: str):
    """
    takes a MatrixTable of annotated Clinvar Variants
    re-annotates these loci with the latest in-house decisions
    reduces dataset to only pathogenic

    re-indexes the data to be queryable on Transcript and Codon
    writes the resulting Table to the specified path

    This was prototyped and executed in a notebook

    Args:
        mt (): MatrixTable of Clinvar variants with VEP anno.
        write_path (): location to write new file to
        new_decisions (str): path to a local clinvar re-summary file

    Returns:
        N/A - write completes within method
    """

    # 1. re-annotate to use the latest decisions
    logging.info(f'loading private clinvar annotations from {new_decisions}')
    clinvar_ht: hl.Table = hl.read_table(new_decisions)

    # remove all non-pathogenic sites
    clinvar_ht = clinvar_ht.filter(
        clinvar_ht.clinvar_significance.lower().contains(PATHOGENIC)
    )

    # minimise MT content
    mt = mt.rows()
    mt = mt.select(tx_csq=mt.vep.transcript_consequences)

    # 1. reduce the annotated clinvar to a minimal representation
    clinvar_ht = mt.rows().join(clinvar_ht)

    # 1. split rows out to separate consequences
    clinvar_ht = clinvar_ht.explode_rows(clinvar_ht.tx_csq)

    # 2. filter down to rows with the relevant content, as a Table
    # i.e. remove non-protein changes, and non-clinvar
    clinvar_ht = clinvar_ht.filter(
        clinvar_ht.tx_csq.consequence_terms.contains('missense_variant')
    )

    # 3. squash the clinvar and protein content into strings
    clinvar_ht = clinvar_ht.annotate(
        clinvar_entry=hl.str('::').join(
            [
                hl.str(clinvar_ht.allele_id),
                clinvar_ht.clinical_significance,
                hl.str(clinvar_ht.gold_stars),
            ]
        ),
        newkey=hl.str('::').join(
            [
                clinvar_ht.tx_csq.protein_id,
                hl.str(clinvar_ht.tx_csq.protein_start),
            ]
        ),
    )

    # create a new table keyed on transcript&position
    temp = clinvar_ht.key_by(clinvar_ht.newkey)
    temp = temp.select(temp.clinvar_entry).collect_by_key()
    temp = temp.annotate(values=hl.set(hl.map(lambda x: x.clinvar_entry, temp.values)))
    temp.write(write_path, overwrite=True)
