#!/usr/bin/env bash

set -ex

# set the date, or provide a default
DATE=${1:-$(date +%F)}
# make a randomized config name
CONFIG_PATH=hail-az://sevgen002sa/test/config-$(LC_ALL=C tr -dc A-Za-z0-9 </dev/urandom | head -c 8).toml

python3 generate_workflow_config.py \
  --dataset severalgenomes \
  --access_level test \
  --driver_image azcpg001acr.azurecr.io/cpg-common/images/cpg_workflows:latest \
  --output_prefix "reanalysis/${DATE}" \
  --extra_datasets severalgenomes rgp \
  --extra_configs reanalysis_global.toml reanalysis_cohort.toml \
  -o ${CONFIG_PATH}

export CPG_CONFIG_PATH=${CONFIG_PATH}
python3 interpretation_runner.py \
  -i hail-az://sevgen002sa/test/reanalysis/986d792a448c66a8a5cfba65434e7d1ce9b1ff_1051-validation.mt \
  --pedigree hail-az://sevgen002sa/test/reanalysis/pedigree.fam \
  --skip_annotation
