#!/usr/bin/env bash

set -ex

COHORT=$1
SUBSAMPLE=${2:-"80"}
TODAY=$(date +%F)

# generate the folder
mkdir -p "inputs/$COHORT"

# generate the pedigree
python helpers/pedigree_from_sample_metadata.py \
    --project "${COHORT}" \
    --output "inputs/${COHORT}/${TODAY}_${SUBSAMPLE}_percent" \
    --plink --hash_threshold \
    "${SUBSAMPLE}"
