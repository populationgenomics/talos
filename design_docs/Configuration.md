# Configuration

The configuration file is a segmented JSON file, containing groups of key: value pairs.

Each of the labelled groups has a `description` label, giving the context for the contained values. The current config groups are described below.

## `moi_tests`

This section contains relates to the application of Mode Of Inheritance tests.

* gnomad_dominant - the max gnomad population AF permitted for Mono-allelic disease
* gnomad_max_homs_dominant - max number of Homs in Gnomad for Mono-allelic disease
* gnomad_max_homs_recessive - max number of Homs in Gnomad for Bi-allelic disease
* gnomad_max_homs_recessive - max number of variants in Gnomad for Mono-allelic disease
* category_2_new_only - analysis flag, if True the category 2 check will only consider where a gene is new in the PanelApp data (i.e. the more granular check of green genes with update MOI is not done)

## `variant_object`

This section relates to the creation of variant objects from CyVCF2 classes.

* csq_string - this is used in Hail to construct the CSQ string, in BCFTools to write into the header, and in downstream analysis to digest it
* var_info_keep - values from the VCF variant INFO field which will be kept

## `filter`

This section relates to the variant filtration and categorisation process, currently used in Hail methods

* af_semi_rare - a minimum AF threshold to be applied to all variants
* min_samples_to_ac_filter - minimum number of samples in the joint-call in order to apply common-in-cohort filter
* ac_threshold - where the joint-call filter is applied, remove all variants more common than this in the callset
* ref_genome - ref_genome
* useless_csq - list of VEP consequences filtered out prior to classification steps
* critical_csq - list of VEP consequences treated as high priority/high impact
* in_silico - separate sub-group, thresholds to apply for each *in silico* tool

## `output`

This section relates to the rubbish HTML document created as output

* seqr_lookup - location of a cache file containing the CPG_ID: SEQR_family_ID mapping
* colours - Hex codes for text colours in the final HTML document
