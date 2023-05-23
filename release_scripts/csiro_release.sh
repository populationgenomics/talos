#!/usr/bin/env bash

# copy ClinVar Data
gsutil -m -u acute-care-321904 cp -r \
    gs://cpg-acute-care-test/csiro
    gs://cpg-acute-care-release/csiro_clinvar_23_05_2023
