# Large Files

For the Stub workflows, the Nextflow configuration expects to find a few files in this folder.

## Annotation Workflow

1. A reference genome matching your input data, in FASTA format, e.g. from `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.{fasta,fasta.fai,dict}`.
2. An Echtvar reference file from https://zenodo.org/records/15222100. Rename this to match the `params.gnomad_zip` entry in the [annotation.config](nextflow/annotation.config) file. This file does not need to be unzipped!
3. An Ensembl GFF3 file, e.g. `Homo_sapiens.GRCh38.113.gff3.gz` from the [Ensembl FTP site](https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens).
4. A MANE Summary text file, e.g. `MANE.GRCh38.v1.4.summary.txt.gz` from the [RefSeq MANE FTP site](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4).
5. Predictions for AlphaMissense, e.g. `AlphaMissense_hg38.tsv.gz` from https://zenodo.org/records/8208688.

The Echtvar reference file in particular has been generated through a compute-intensive encoding of select gnomAD V4.1
raw data. If you have a preferred tool internally for applying gnomAD annotations to a VCF, you can use that instead,
provided it can add the following fields to the VCF `INFO` column:

- `gnomad_AC_joint` (`int`)
- `gnomad_AF_joint` (`float`)
- `gnomad_AC_joint_XY` (`int`)
- `gnomad_HomAlt_joint` (`int`)

## Talos Workflow

1. `genes_to_phenotypes.txt` from https://hpo.jax.org/data/annotations
2. `hp.obo` from https://hpo.jax.org/data/ontology
3. `phenio.db.gz` from https://data.monarchinitiative.org/monarch-kg/latest/phenio.db.gz
4. `clinvarbitration` data from Zenodo https://zenodo.org/records/15896156

> **NOTE** the ClinvArbitration data is updated and re-uploaded monthly. A Stub has been provided in this repository, matching the dummy test data (`nextflow/inputs/clinvarbitration.tar.gz`) but that is not suitable for real analyses (only contains 4 variants).
