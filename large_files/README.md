# Large Files

Talos requires a number of large files, for both the minimal test workflow, and real-world analyses. These files are too large to store in GitHub, so they are not included in this repository. Instead, a file [gather_files.sh](gather_files.sh) is provided, which will download the required files from various sources prior to running analyses.

The `gather_files.sh` script will download the files to the local directory, wherever it is run. That can either be this `large_files` directory, or any other directory you choose. When running the [Annotation](../nextflow/annotation.nf) and [Talos](../nextflow/talos.nf) workflows, you can override the default `large_files` directory by setting the `--large_files` parameter when invoking the nextflow run command, indicating the directory where you have downloaded the files.

The file names created by the `gather_files.sh` script are the same as the default names in the [annotation.config](../nextflow/annotation.config) and [talos.config](../nextflow/talos.config) files - if file names in the configuration files are altered, the file names should be udpated the match.

> **n.b.** Talos has only been tested on GRCh38, and in that particular reference genome is used in a number of places. The download script fetches a copy of HG38 as part of the setup process, and Talos is not technically supported on other reference genomes, so your variant data will have to be aligned to GRCh38 for the time being. If you have a need to run Talos on a different reference genome, please raise an issue in the GitHub repository.

## Annotation Workflow

1. A reference genome matching your input data, in FASTA format, e.g. from `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.{fasta,fasta.fai,dict}`. This is not downloaded by [gather_files.sh](gather_files.sh), as this may be specific to your input data.
2. An Echtvar reference file from https://zenodo.org/records/15222100.
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
