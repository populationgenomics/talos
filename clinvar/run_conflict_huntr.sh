#!/usr/bin/env bash

set -ex

# set the date, or provide a default
DATE=${1:-$(date +%F)}

# run
analysis-runner \
    --config clinvar/config.toml \
    --dataset acute-care \
    --description "Process clinvar annotations" \
    -o "clinvar_${DATE}" \
    --access-level test \
    --cpu 8 \
        python3 clinvar/conflict_huntr.py \
        -s gs://cpg-acute-care-test/clinvar_annotations/submission_summary.txt.gz \
        -v gs://cpg-acute-care-test/clinvar_annotations/variant_summary.txt.gz \
        -o clinvar_reprocessed.ht
