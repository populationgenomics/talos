# Helpers

The [helpers](../helpers) directory is designed to be a collection place for single purpose scripts, useful
in setting up or configuring the analysis, but not a direct part of the end-to-end workflow.

## Run Setup with [prepare_aip_cohort.py](../helpers/prepare_aip_cohort.py)

To simplify setup and running of AIP for CPG cohorts the script has been included. This is
a multifunctional script, and collects a number of behaviours into a single script, at the
expense of being specific to CPG usage.

This script requires usage of a config file, containing a subsection for the relevant cohort. See
[the config example](../reanalysis/reanalysis_global.toml) `cohorts` subsections for details.

1. Queries [Metamist](https://sample-metadata.populationgenomics.org.au/), the CPG metadata server
2. Uses the results to generate a pedigree for the cohort being processed
3. [Combines HPO metadata](HPO_Panel_Matching.md) with a Human Phenotype Ontology file database to
   find any phenotype-matched panels to use for each family
4. [optional] Ingests Seqr metadata for the same cohort, used to match sample and family IDs between
   the analysis files (internal IDs) and Seqr (auto-inc. upon ingestion)
5. [optional] Takes a flag to generate info as an `exomes` cohort

The output of this script creates a number of outputs:

1. A folder `inputs/<cohort>` containing output files
2. A pedigree file `pedigree.ped`
3. A mapping of internal to external IDs `external_lookup.json`
4. A mapping of participants to phenotype-matched panels `participant_panels.json`
5. A processed Seqr family lookup `seqr_processed.json`
6. A cohort TOML file containing settings to use `cohort_config.toml`
7. The pedigree, lookup, seqr mapping, and panel data all copied to the relevant GCP bucket
8. [optional] If the exome flag is used, this will use `_exomes` in file paths to prevent clashes
   with the genome data

## Central Index Generator [report_hunter.py](../helpers/report_hunter.py)

This script crawls all project GCP buckets it has access to and identifies the latest
HTML reports for each cohort. These collected reports are assembled into an index HTML
page, which is posted to the `common` bucket.

The source of the project list to check is Metamist, using the projects end point. The
buckets are listed using a `google.storage` method, accumulating all report paths in
the relevant buckets, which are then sorted by modified datetime to find latest.

Reports are sub-divided by Genome/Exome, and Familial(default)/Singletons, with the
latest report in each category hyperlinked to the index page. This page is then written
to [popgen.rocks/aip_index](http://popgen.rocks/aip_index).

Due to transient cohort dependencies, any user with ability to read any one of the CPG
projects will have permission to view the index page, though they may not have permission
to view all linked reports.

## Forbidden Gene Finder [forbidden_gene_check.py](../helpers/forbidden_gene_check.py)

This script exists to enable simulated retrospective analysis, by asking the question:

"At a given point in time, which genes were green PanelApp, and which genes were added
to Panels between that date and now."

```commandline
Given a list of phenotype-matched panels, and a number of dates

Find the current state of PanelApp (all green gene symbols on the aggregated panel list)

For each date in turn find the aggregate gene list of all relevant panels at the given
point, and find a diff between that list and the current gene list.

Using the blacklisted/forbidden gene list mechanic, we can manually cut genes out of the
analysis where we would not have seen them, had we queried PanelApp at that point in time.
```
