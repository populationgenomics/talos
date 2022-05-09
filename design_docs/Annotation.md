# Annotation

The code designed by @vladsaveliev to apply annotations using Hail Batch has been
extracted here from the core CPG pipeline.

## Input

The annotation process takes data in a VCF format, partitioning data based on GATK-suggested
intervals to parallelise the process across a number of separate Hail Batch jobs.

Each Job takes a region of the VCF and annotates each in turn using command-line VEP. These
intermediate annotation files are combined into a single Hail Table.

The VCF is then read into MatrixTable format, where the Table data is applied into to the MT.
The product of these stages is the variant data in a


## Reference Material

The annotation process makes use of multiple different resources:

- VEP provides variant consequence annotations, and LOFTEE annotations through a plugin
- ClinVar, Gnomad, ExAC, and a number of in silico tool annotations are provided directly
via annotation from existing Hail Tables

The VEP input data is contained within a single directory. For CPG uses, this is stored in
`gs://cpg-reference/hg38`, which is a copy of `gs://gcp-public-data--broad-references/hg38`.

Hail Table-format annotations are stored in `gs://cpg-seqr-reference-data`, which is a copy
of `gs://seqr-reference-data`; This contains a table for ClinVar annotations, and another
table featuring a collection of annotations from various resources. For CPG use these are
copied from the Broad Institute-prepared data, with the scripts used to generate them stored
[here](https://github.com/broadinstitute/hail-elasticsearch-pipelines/tree/master/download_and_create_reference_datasets/v02).
Specifically, [this is used to generate the combined data table](https://github.com/broadinstitute/hail-elasticsearch-pipelines/blob/master/download_and_create_reference_datasets/v02/hail_scripts/write_combined_reference_data_ht.py),
and [this is used to generate new dumps of clinvar data](https://github.com/broadinstitute/hail-elasticsearch-pipelines/blob/master/hail_scripts/utils/clinvar.py).
