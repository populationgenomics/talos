#!/usr/bin/env bash

set -x

# 'canonical' files
CANONICAL_FILES="gs://cpg-acute-care-test/reanalysis/pre_panelapp_mendeliome.json gs://cpg-acute-care-test/reanalysis/vqsr_header_line.txt gs://cpg-acute-care-test/reanalysis/csq_header_line.txt"
COHORT_LIST=("acute-care")

for cohort in ${COHORT_LIST[@]}; do
    gsutil -m cp ${CANONICAL_FILES} "gs://cpg-${cohort}-test/reanalysis"
done
