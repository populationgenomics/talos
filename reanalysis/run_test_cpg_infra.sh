#!/usr/bin/env bash

set -ex

analysis-runner \
  --dataset severalgenomes \
  --description "MS CPG infra test" \
  --access-level test \
  --output-dir "miah_test"
  reanalysis/test_cpg_infra.py
    --blob hail-az://sevgen002sa/cpg-severalgenomes-test/hello.txt