#!/usr/bin/env bash
# https://github.com/populationgenomics/seqr-private/issues/11

# copy single files
gsutil -m -u acute-care-321904 cp \
    gs://cpg-acute-care-test/csiro_copy/callset.vcf.bgz \
    gs://cpg-acute-care-test/csiro_copy/hpo_terms.json \
    gs://cpg-acute-care-test/csiro_copy/hpo_terms.tsv \
    gs://cpg-acute-care-test/csiro_copy/pedigree.fam \
    gs://cpg-acute-care-release/csiro
