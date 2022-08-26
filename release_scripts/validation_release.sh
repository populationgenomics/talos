#!/usr/bin/env bash
# https://github.com/populationgenomics/automated-interpretation-pipeline/issues/101

# copy single files
gsutil -m -u validation-351902 cp \
    gs://cpg-validation-test/reanalysis/vqsr_header_line.txt \
    gs://cpg-validation-test/reanalysis/pre_panelapp_mendeliome.json \
    gs://cpg-validation-test/reanalysis/csq_header_line.txt \
    gs://cpg-validation-test/reanalysis/pedigree.fam \
    gs://cpg-validation-test/reanalysis/reanalysis_conf.json \
    gs://cpg-validation-test/reanalysis/external_lookup.json \
    gs://cpg-validation-release/reanalysis

# recursive copy of MT
gsutil -m -u validation-351902 cp -r \
    gs://cpg-validation-main/mt/986d792a448c66a8a5cfba65434e7d1ce9b1ff_1051-validation.mt \
    gs://cpg-validation-release/reanalysis/2022-09-19_validation.mt

# copy and rename
gsutil -m -u validation-351902 cp \
    gs://cpg-validation-test/reanalysis/2022-08-19/summary_results.json \
    gs://cpg-validation-release/reanalysis/2022-09-19_expected.json
