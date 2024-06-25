# Configuration File

The configuration file is a segmented TOML file, containing groups of key: value pairs. The current config groups are
described below. The strategy for digesting the config file is based on the core CPG's config management approach,
explained [here](https://github.com/populationgenomics/team-docs/blob/main/cpg_utils_config.md)

## `workflow`

In all CPG config files, `workflow` is the first section. This section contains the core settings for the workflow, and
is used to store details of the current setup, e.g. the authenticated user/project:

* presentation - this is used to select the presentation mode for the current run. Only the CPG mode currently exists,
  which will trigger the creation of a HTML summary based on our templates
* status_reporter - the status reporter to choose, mirroring CPG's core pipeline default. `status reporting` relates to
  how and where results of a completed run are logged. `metamist` is the only implemented reporter in this repository,
  which will write analysis entries to the CPG Metamist metadata database
* ignore_categories - a list, add category titles here to remove them from consideration

## `panels`

The `panels` section relates to the gene panel querying phase

* panelapp - the URL of the panelapp instance to use in this analysis
* default_panel - the base panel for this analysis (defaults to the Mendeliome, 137)
* require_pheno_match - genes to remove from the base panel (for noise reasons), but to permit if they are present on a
  phenotype-matched panel

## clinvar

This section relates to the ClinVar re-summary and annotation process. Small.

* filter_benign - a list of submission sites where we're happy to disregard benign submissions, e.g. population
  databases with insufficient per-sample context

## `moi_tests`

This section contains relates to the application of Mode Of Inheritance tests. These are not set based on the prevalence
any individual disease in any specific population. Feedback from clinical teams will be used to further refine these
settings

* gnomad_dominant - the max gnomad population AF permitted for Mono-allelic disease
* gnomad_max_homs_dominant - max number of Homs in Gnomad for Mono-allelic disease
* gnomad_max_homs_recessive - max number of Homs in Gnomad for Bi-allelic disease
* gnomad_max_ac_dominant - max number of variants in Gnomad for Mono-allelic disease
* category_2_new_only - analysis flag, if True the category 2 check will only consider where a gene is new in the
  PanelApp data (i.e. the more granular check of green genes with update MOI is not done)

## `csq`

This section determines which VEP elements are retained in the VCF, and are used in constructing the intermediate VCF
header, and decoding the VCF data.

* csq_string - array of strings, describing the content in the CSQ field

## `filter`

This section relates to the variant filtration and categorisation process, currently executed in Hail query. These
default settings are primarily based on a clinical-team provided specification, and can be replaced for a bespoke run.

* af_semi_rare - a minimum AF threshold to be applied to all variants
* min_samples_to_ac_filter - minimum number of samples in the joint-call in order to apply common-in-cohort filter
* ac_threshold - where the joint-call filter is applied, remove all variants where the joint-call AC/AD (variant
  calls/total alleles) is above this value (default cut-off is AC >= 10%)
* ref_genome - ref_genome
* useless_csq - list of VEP consequences filtered out prior to classification steps
* critical_csq - list of VEP consequences treated as high priority/high impact
* in_silico - separate sub-group, thresholds to apply for each *in silico* tool. These parameters are taken
  from [work by *Pejaver et al.*](https://www.biorxiv.org/content/10.1101/2022.03.17.484479v1)

## `categories`

This section contains the categories included in the current codebase, and a brief description of each. As new
categories are added, this should be kept updated.

## `output`

This section relates to the HTML document created as output

* seqr_lookup - location of a cache file containing the CPG_ID: SEQR_family_ID mapping
* colors - Hex codes for text colours in the final HTML document

## `cohorts`

A section containing all the cohorts known about, and any data specific to each

* solved_cases - a list of any Proband/Family IDs which we want to remove from the final report
* cohort_panel

## `cohorts.sequencing_type`

This is a sub-section, containing data relevant to a specific sequencing type within a cohort.

* external_labels - path to a file containing labels to apply to variants in individuals (e.g. failed validation)
* seqr_project - the project ID for this cohort in SEQR
* seqr_lookup - the mapping of CPG internal ID to Seqr ID
* seqr_instance - URL of Seqr instance containing this cohort
