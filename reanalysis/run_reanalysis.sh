#!/usr/bin/env bash

set -ex

# set the date, or provide a default
PAP_DATE=${1:-"2021-09-03"}

# run
analysis-runner \
  --dataset acute-care \
  --description "run reanalysis draft" \
  -o "reanalysis/${PAP_DATE}" \
  --access-level test \
  reanalysis/wrapper.py \
    --config_json gs://cpg-acute-care-test/reanalysis/reanalysis_conf.json \
    --matrix gs://cpg-acute-care-main/mt/acute-care.mt \
    --panel_genes gs://cpg-acute-care-test/reanalysis/pre_panelapp_mendeliome.json
