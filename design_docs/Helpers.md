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

### Pedigree manipulation

The Pedigree data pulled directly from the SM-API contains all members and relationships in the context
of the participants' external identifiers (PID). This can't be matched up with the VCF contents, as in
sequence data files processing is instead done in the context of the Sample IDs (SID).

To complicate things, there is no guarantee of a 1:1 relationship between participants and samples (e.g.
replacement samples, but also Genome and Exome library preparations would be considered separate sequences
derived from the same physical sample).
To overcome any potential ambiguity, this process uses the following steps:

- obtain the Pedigree data with external identifiers
- obtain the mapping of PID:SID, allowing for multiple SID/PID
- dump a dictionary to file, allowing lookup of all SID -> PID
- for each PID, collect all possible SID as a list
- when writing the translated pedigree:
  - use itertools.product to generate intra-family combinations of the proband, father, & mother e.g.
  - 1 sample per participant == 1^3 = 1 row
  - 2 samples per participant == 2^3 = 8 rows for all possible combinations

This pedigree will encompass all possible sample collections, and with a high level of sample replacement
could expand massively. For an individual use case sample IDs could be taken from a joint-call or sample-list,
and used to reduce the number of rows read from this PED to only those relevant to a current analysis
