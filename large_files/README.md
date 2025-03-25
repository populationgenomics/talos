## Large Files

For the Stub workflows, the Nextflow configuration expects to find a few files in this folder:

1. A reference genome matching your input data, in FASTA format
2. An echtvar reference file - we have a pre-generated file we are able to share, please get in touch
3. `genes_to_phenotypes.txt` from https://hpo.jax.org/data/annotations
4. `hp.obo` from https://hpo.jax.org/data/ontology
5. `phenio.db.gz` from https://data.monarchinitiative.org/monarch-kg/latest/phenio.db.gz

The Echtvar reference file in particular has been generated through a compute-intensive encoding of select gnomAD V4.1
raw data. If you have a preferred tool internally for applying gnomAD annotations to a VCF, you can use that instead,
provided it can add the following fields to the VCF INFO column:
    - gnomad_AC_joint (int)
    - gnomad_AF_joint (float)
    - gnomad_AC_joint_XY (int)
    - gnomad_HomAlt_joint (int)
