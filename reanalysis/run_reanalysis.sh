#!/usr/bin/env bash

set -ex

# set the date, or provide a default
DATE=${1:-$(date +%F)}

analysis-runner \
  --config reanalysis/reanalysis_global.toml \
  --config reanalysis/reanalysis_cohort.toml \
  --image azcpg001acr.azurecr.io/cpg-common/images/cpg_aip \
  --dataset rgp \
  --description "AIP run" \
  -o "reanalysis_train/${DATE}" \
  --access-level test \
  reanalysis/interpretation_runner.py \
    -i hail-az://raregen001sa/test/inputs/rgp/rgp_train.vcf.bgz \
    --pedigree hail-az://raregen001sa/test/inputs/rgp/rgp_train.fam \
    --participant_panels hail-az://raregen001sa/test/inputs/rgp/rgp_default_panels.json
