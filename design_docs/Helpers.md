# Helpers

The [helpers](../helpers) directory is designed to be a collection place for single purpose scripts, useful
in setting up or configuring the analysis, but not a direct part of the end-to-end workflow.

## Pedigree

Running accurate genetic analysis is heavily dependent on an accurate representation of relationships between samples.
For the CPG, this data is typically stored within the sample-metadata database, accessible through the corresponding API
client.

The [Pedigree-generation](../helpers/pedigree_from_sample_metadata.py) script has been included which
will query the sample-metadata API for a given cohort, and rearrange the results into two outputs:

- A Pedigree
- A mapping of the external individual IDs to CPG sample IDs

Additional options can be supplied to:

- export the Pedigree as PLINK format instead of .ped format
- export the Pedigree with all individuals as singletons
