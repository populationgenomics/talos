#!/usr/bin/env bash

set -ex

# set the date, or provide a default
DATE=${1:-$(date +%F)}

analysis-runner \
  --config reanalysis/reanalysis_global.toml \
  --config reanalysis/reanalysis_cohort.toml \
  --image azcpg001acr.azurecr.io/cpg-common/images/cpg_workflows \
  --dataset severalgenomes \
  --description "AIP run" \
  -o "reanalysis/${DATE}" \
  --access-level test \
  reanalysis/interpretation_runner.py \
    -i hail-az://sevgen002sa/test/reanalysis/2022-11-16/prior_to_annotation.vcf.bgz \
    --pedigree hail-az://sevgen002sa/test/reanalysis/pedigree.fam