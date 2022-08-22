#!/usr/bin/env bash

set -ex

# set the date, or provide a default
PAP_DATE=${1:-$(date +%F)}

# run
analysis-runner \
  --dataset acute-care \
  --description "Run Comparison" \
  -o "reanalysis/comparison/${PAP_DATE}" \
  --access-level test \
  comparison/comparison_wrapper.py \
    --results_folder gs://cpg-acute-care-test/reanalysis/2022-08-19 \
    --seqr gs://cpg-acute-care-test/reanalysis/comparison/seqr_acute_care_tags.tsv \
    --mt gs://cpg-acute-care-main/mt/986d792a448c66a8a5cfba65434e7d1ce9b1ff_1051-acute-care.mt
