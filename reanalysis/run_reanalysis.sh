#!/usr/bin/env bash

set -ex

# set the date, or provide a default
DATE=${1:-$(date +%F)}

analysis-runner \
  --config reanalysis/reanalysis_global.toml \
  --config reanalysis/reanalysis_cohort.toml \
  --dataset acute-care \
  --description "AIP run" \
  -o "reanalysis/${DATE}" \
  --access-level test \
  reanalysis/interpretation_runner.py \
    --input_path gs://cpg-acute-care-test/reanalysis/2011-11-11/prior_to_annotation.vcf.bgz \
    --panel_genes gs://cpg-acute-care-test/reanalysis/pre_panelapp_mendeliome.json \
    --plink_file gs://cpg-acute-care-test/reanalysis/acute-care-plink.fam
