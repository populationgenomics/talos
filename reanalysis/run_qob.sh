#!/usr/bin/env bash

set -ex

export CPG_DEPLOY_CONFIG
analysis-runner \
  --dataset severalgenomes \
  --description "MS QoB test" \
  --access-level test \
  --output-dir "miah_test" \
  --env CPG_CONFIG_PATH="hail-az://sevgen002sa/cpg-severalgenomes-main/cpg-config.toml" \
  reanalysis/test_qob.py