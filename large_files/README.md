# Large Files

Talos requires a number of large files, for both the minimal test workflow, and real-world analyses. These files are too large to store in GitHub, so they are not included in this repository. Instead, a file [gather_files.sh](gather_files.sh) is provided, which will download the required files from various sources prior to running analyses.

The `gather_files.sh` script will download the files to the local directory, wherever it is run. That can either be this `large_files` directory, or any other directory you choose. When running the [Annotation](../nextflow/annotation.nf) and [Talos](../nextflow/talos.nf) workflows, you can override the default `large_files` directory by setting the `--large_files` parameter when invoking the nextflow run command, indicating the directory where you have downloaded the files.

The file names created by the `gather_files.sh` script are the same as the default names in the [annotation.config](../nextflow/annotation.config) and [talos.config](../nextflow/talos.config) files - if file names in the configuration files are altered, the file names should be udpated to match. The exception is the reference genome

If you already have any of these data files or reference genomes present locally, e.g. from a prior deployment or other local work, the script will not re-download them, so you can run the script multiple times without worrying about overwriting existing files. This behaviour requires the file names to match exactly (e.g. the reference genome is expected as `ref.fa`).

> **n.b.** Talos has only been tested on GRCh38, and in that particular reference genome is used in a number of places. The download script fetches a gzip-compressed copy of HG38 as part of the setup process, and Talos is not technically supported on other reference genomes, so your variant data will have to be aligned to GRCh38 for the time being. If you have a need to run Talos on a different reference genome, please raise an issue in the GitHub repository. To use the reference genome downloaded by this scripted process, you will need to decompress with `gunzip` and index with `samtools faidx` prior to its use in the workflow.

## Files

The `gather_files.sh` script downloads the following:

1. A compressed GRCh38 reference genome, in gzipped FASTA format, e.g. from `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/` (see note above about decompressing and indexing)
2. An Echtvar reference file from `https://zenodo.org/records/15222100`
3. An Ensembl GFF3 file, e.g. `Homo_sapiens.GRCh38.113.gff3.gz` from the [Ensembl FTP site](https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens).
4. A MANE Summary text file, e.g. `MANE.GRCh38.v1.4.summary.txt.gz` from the [RefSeq MANE FTP site](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4).
5. Predictions for AlphaMissense, e.g. `AlphaMissense_hg38.tsv.gz` from https://zenodo.org/records/8208688.
6. `genes_to_phenotypes.txt` and `hp.obo` from `https://github.com/obophenotype/human-phenotype-ontology/releases`
7. `phenio.db` from https://data.monarchinitiative.org/monarch-kg/latest/phenio.db.gz, then decompressed
8. `ClinvArbitration` data from Zenodo https://zenodo.org/records/16792026

> **NOTE** the ClinvArbitration data is updated and re-uploaded monthly. This should be downloaded from Zenodo prior to attempting to run Talos, and an updated version should be downloaded each month to stay current.
