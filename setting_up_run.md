# setting up a new run

CPG-specific guidance on how to generate input files for a new run

## Initial checks

1. Check that
   the [token map](https://github.com/populationgenomics/cpg-infrastructure-private/blob/main/tokens/repository-map.json)
   allows AIP to run against the relevant project's bucket & contents
2. Check that you have permission to query Metamist for the relevant project (e.g. you can run queries for the project
   on the [Swagger page](https://sample-metadata.populationgenomics.org.au/swagger))
3. Check that you have permission to upload relevant files to the project bucket (current process involves uploading
   files to the permissive `-test` bucket for ease)

## Setup Shortcut

A Helper script has been written called [prepare_aip_cohort.py](helpers/prepare_aip_cohort.py), which will generate all
the required inputs for the analysis run, and copy them directly to GCP. To use this script:

1. To add Seqr data automatically, you will need to obtain the seqr metadata file as described below, and save it to a
   JSON file.
2. Run the script in the following way:

* set the environment variable `CPG_CONFIG_PATH`, pointing to a local config,
  e.g. [reanalysis/reanalysis_global.toml](reanalysis/reanalysis_global.toml)
* set the argument `--project` to match the metamist cohort name
* `--seqr` should point to the seqr metadata file, if required
* `--obo` should point to the HPO obo file, if required
* `-e` is required if the samples required are exomes
* `--singletons` is required if the samples should be processed as singletons

e.g.

```bash
CPG_CONFIG_PATH=reanalysis/reanalysis_global.toml \
    python helpers/prepare_aip_cohort.py \
    --project broad-rgp \
    --seqr inputs/broad-rgp/seqr.json \
    --obo helpers/hpo_terms.obo
```

This script should generate all the required files, and copy some to the relevant GCP bucket in `<root>/reanalysis`. The
generated files are:

* Pedigree file representing this cohort
* Processed Seqr data (if provided)
* Mapping of CPG-to-External IDs
* The per-participant panels to use
* Pre-PanelApp Mendeliome file (if relevant)
* Cohort-specific config to use (containing above filepaths)

## Seqr Family Lookup

* Optional, if provided the final HTML will hyperlink to Seqr variants and families
* Go to a relevant project page in Seqr (clicking the project name)
* (in Chrome) Burger -> More Tools -> Developer Tools (Opt+Cmd+I) & Open the Network tab
* Refresh, then find the named request `get_overview`
* On the right-side of the console, find the `Response` tab, which should show some JSON
* Select all this JSON and save to a file

## Starting the run

Create a run script based on the example [run_reanalysis.sh](reanalysis/run_reanalysis.sh)

All paths referenced here must be GCP bucket paths, as no files will be copied to GCP at runtime. For input genomic data
this will involve finding the appropriate paths to use. For input files this will involve transcribing the paths you
previously created into this script as inputs.
