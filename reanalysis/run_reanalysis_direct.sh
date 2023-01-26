#!/usr/bin/env bash

set -ex

# Make sure to export CPG_DEPLOY_CONFIG first?

# set the date, or provide a default
#DATE=${1:-$(date +%F)}
DATE="2023-01-26"
# make a randomized config name
CONFIG_PATH=hail-az://sevgen002sa/test/config-$(LC_ALL=C tr -dc A-Za-z0-9 </dev/urandom | head -c 8).toml

  # --deploy_config ~/sources/cpg/cpg-deploy/azure/deploy-config.prod.json \
  # --server_config ~/sources/cpg/cpg-deploy/aip/terraform.tfvars.json \
python3 reanalysis/generate_workflow_config.py \
  --dataset rgp \
  --access_level test \
  --driver_image azcpg001acr.azurecr.io/cpg-common/images/cpg_aip:latest \
  --output_prefix "reanalysis_train/${DATE}" \
  --extra_datasets severalgenomes rgp \
  --image_base azcpg001acr.azurecr.io/cpg-common/images \
  --reference_base hail-az://azcpg001sa/reference \
  --extra_configs reanalysis/reanalysis_global.toml reanalysis/reanalysis_cohort.toml \
  -o ${CONFIG_PATH}
  
export CPG_CONFIG_PATH=${CONFIG_PATH}
python3 reanalysis/interpretation_runner.py \
  -i hail-az://raregen001sa/test/inputs/rgp/rgp_train.vcf.bgz \
  --pedigree hail-az://raregen001sa/test/inputs/rgp/rgp_train.fam \
  --participant_panels hail-az://raregen001sa/test/inputs/rgp/rgp_default_panels.json
