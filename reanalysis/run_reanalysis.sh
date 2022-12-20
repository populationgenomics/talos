#!/usr/bin/env bash

set -ex

# set the date, or provide a default
DATE=${1:-$(date +%F)}

analysis-runner \
  --config reanalysis/reanalysis_global.toml \
  --config reanalysis/reanalysis_cohort.toml \
  --image azcpg001acr.azurecr.io/cpg-common/images/cpg_workflows \
  --dataset rgp \
  --description "AIP run" \
  -o "reanalysis/${DATE}" \
  --access-level test \
  reanalysis/interpretation_runner.py \
    -i hail-az://raregen001sa/test/inputs/joint-called-vcf_20221114/RGP_subset_samples.vcf.bgz \
    --pedigree hail-az://raregen001sa/test/inputs/joint-called-vcf_20221114/RGP_Cases_for_MSFT_AIP_v0_trial.xlsx.fam