#!/usr/bin/env bash

set -ex

# set the date, or provide a default
PAP_DATE=${1:-"2021-09-03"}

# run
analysis-runner \
  --dataset acute-care \
  --description "run automated_interpretation_pipeline process" \
  -o "reanalysis/${PAP_DATE}" \
  --access-level test \
  reanalysis/interpretation_runner.py \
    --config_json gs://cpg-acute-care-test/reanalysis/reanalysis_conf.json \
    --input_path gs://cpg-acute-care-test/reanalysis/${PAP_DATE}/hail_105_ac.mt/ \
    --panel_genes gs://cpg-acute-care-test/reanalysis/pre_panelapp_mendeliome.json \
    --pedigree gs://cpg-acute-care-test/reanalysis/acute-care-plink.fam
