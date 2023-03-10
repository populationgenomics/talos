"""
method file for re-sorting clinvar annotations by codon
utilises the latest in-house clinvar summary data
this should probably be transformed into a broader wrapper;

1. download ClinVar VCF
2. annotate that VCF with VEP consequences
3. substitute out the clinvar decisions with CPG versions
4. munge the data to index all CPG decisions on protein change
5. write resulting table out

The main reason not to do this just yet is that we won't need
to re-run it frequently, and hopefully it won't be long before
VEP-on-Hail-Query works, in which case that whole workflow
would be a trivial extension from this
"""


import logging

import hail as hl

from reanalysis.hail_filter_and_label import PATHOGENIC, ONE_INT


def protein_indexed_clinvar(mt: hl.MatrixTable, write_path: str, new_decisions: str):
    """
    takes a MatrixTable of annotated Clinvar Variants
    re-annotates these loci with the latest in-house decisions
    reduces dataset to only pathogenic

    re-indexes the data to be queryable on Transcript and Codon
    writes the resulting Table to the specified path

    Prototyped and executed in a notebook

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

    # minimise clinvar vcf content prior to the join
    mt = mt.rows()
    mt = mt.select(tx_csq=mt.vep.transcript_consequences)

    # join the new decisions with the clinvar allele consequence annotations
    clinvar_ht = mt.rows().join(clinvar_ht)

    # 2. split rows out to separate transcript consequences
    clinvar_ht = clinvar_ht.explode_rows(clinvar_ht.tx_csq)

    # 3. filter down to rows with the relevant content
    # filter for missense & single substitution (prevent frameshift)
    # a reasonable filter here would also include MANE transcripts
    clinvar_ht = clinvar_ht.filter(
        (hl.len(clinvar_ht.alleles[0]) == ONE_INT)
        & (hl.len(clinvar_ht.alleles[1]) == ONE_INT)
        & clinvar_ht.tx_csq.consequence_terms.contains('missense_variant')
    )

    # 4. squash the clinvar and protein content into single strings
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

    # 5. re-key table on transcript & residue
    clinvar_ht = clinvar_ht.key_by(clinvar_ht.newkey)

    # 6. collect all ClinVar annotations at each residue
    clinvar_ht = clinvar_ht.select(clinvar_ht.clinvar_entry).collect_by_key()

    # 7. squash the multiple clinvar entries back to a single string
    clinvar_ht = clinvar_ht.transmute(
        clinvar_alleles=hl.str('+').join(
            hl.set(hl.map(lambda x: x.clinvar_entry, clinvar_ht.values))
        )
    )

    # 8. write the table of all ENSP:residue#: Clinvar[+Clinvar,]
    clinvar_ht.write(write_path, overwrite=True)
