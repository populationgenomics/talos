# Configuration File

The configuration file is a segmented JSON file, containing groups of key: value pairs.

Each of the labelled groups has a `description` label, giving the context for the contained values. The current config groups are described below.

## `moi_tests`

This section contains relates to the application of Mode Of Inheritance tests. These are not set based on the prevalence
any individual disease in any specific population. Feedback from clinical teams will be used to further refine these
settings

* gnomad_dominant - the max gnomad population AF permitted for Mono-allelic disease
* gnomad_max_homs_dominant - max number of Homs in Gnomad for Mono-allelic disease
* gnomad_max_homs_recessive - max number of Homs in Gnomad for Bi-allelic disease
* gnomad_max_ac_dominant - max number of variants in Gnomad for Mono-allelic disease
* category_2_new_only - analysis flag, if True the category 2 check will only consider where a gene is new in the PanelApp data (i.e. the more granular check of green genes with update MOI is not done)

## `variant_object`

This section relates to the creation of variant objects from CyVCF2 classes. When parsing a VCF using Cyvcf2, the
variant representation objects use a custom representation of the INFO fields, which cannot be pickled. This limits the
possible parallelisation of analysis downstream. Extracting out into dictionaries opens up opportunities downstream, and
these fields are used to set up the dictionary fields.

CSQ_STRING is used in 3 separate places:

1. Hail to construct the CSQ string
2. BCFTools to write into the header
3. in downstream analysis to digest the CSQ content

* csq_string - |-delimited series of values, describing the content in the CSQ field

## `filter`

This section relates to the variant filtration and categorisation process, currently executed in Hail query. These
default settings are primarily based on a clinical-team provided specification, and can be replaced for a bespoke run.

* af_semi_rare - a minimum AF threshold to be applied to all variants
* min_samples_to_ac_filter - minimum number of samples in the joint-call in order to apply common-in-cohort filter
* ac_threshold - where the joint-call filter is applied, remove all variants where the joint-call AC/AD (variant calls/total alleles) is above this value (default cut-off is AC >= 10%)
* ref_genome - ref_genome
* useless_csq - list of VEP consequences filtered out prior to classification steps
* critical_csq - list of VEP consequences treated as high priority/high impact
* in_silico - separate sub-group, thresholds to apply for each *in silico* tool. These parameters are taken from [work by *Pejaver et al.*](https://www.biorxiv.org/content/10.1101/2022.03.17.484479v1)

## `output`

This section relates to the rubbish HTML document created as output. To be discontinued ASAP

* seqr_lookup - location of a cache file containing the CPG_ID: SEQR_family_ID mapping
* colors - Hex codes for text colours in the final HTML document
