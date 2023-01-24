#!/usr/bin/env bash

set -ex

# Make sure to export CPG_DEPLOY_CONFIG first?

# set the date, or provide a default
#DATE=${1:-$(date +%F)}
DATE="2023-01-23"
# make a randomized config name
CONFIG_PATH=hail-az://sevgen002sa/test/config-$(LC_ALL=C tr -dc A-Za-z0-9 </dev/urandom | head -c 8).toml

  # --deploy_config ~/sources/cpg/cpg-deploy/azure/deploy-config.prod.json \
  # --server_config ~/sources/cpg/cpg-deploy/aip/terraform.tfvars.json \
python3 generate_workflow_config.py \
  --dataset rgp \
  --access_level test \
  --driver_image azcpg001acr.azurecr.io/cpg-common/images/cpg_aip:latest \
  --output_prefix "reanalysis/${DATE}" \
  --extra_datasets severalgenomes rgp \
  --image_base azcpg001acr.azurecr.io/cpg-common/images \
  --reference_base hail-az://azcpg001sa/reference \
  --extra_configs reanalysis_global.toml reanalysis_cohort.toml \
  -o ${CONFIG_PATH}
  
export CPG_CONFIG_PATH=${CONFIG_PATH}
python3 interpretation_runner.py \
  -i hail-az://raregen001sa/test/inputs/rgp/rgp_test.vcf.bgz \
  --pedigree hail-az://raregen001sa/test/inputs/joint-called-vcf_20221114/RGP_Cases_for_MSFT_AIP_v0_trial.xlsx.fam \
  --participant_panels hail-az://raregen001sa/test/inputs/rgp/rgp_party_panels.json
