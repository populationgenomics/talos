#!/usr/bin/env bash

set -ex

# run
analysis-runner \
  --dataset acute-care \
  --description "run automated_interpretation_pipeline process" \
  -o reanalysis/comp_het_test \
  --access-level test \
  reanalysis/isolated_runner.py \
    --matrix_path gs://cpg-acute-care-test/reanalysis/2021-09-03/hail_categorised.vcf.bgz.mt
