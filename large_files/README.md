# Large Files

Talos requires a number of large files, for both the minimal test workflow, and real-world analyses. These files are too large to store in GitHub, so they are not included in this repository. Instead, a file [gather_files.sh](gather_files.sh) is provided, which will download the required files from various sources prior to running analyses.

The `gather_files.sh` script will download the files to the local directory, wherever it is run. That can either be this `large_files` directory, or any other directory you choose. When running the workflow, you can override the default `large_files` directory by setting `params.large_files` or `--large_files`, indicating the directory containing the files.

The file names created by the `gather_files.sh` script are the same as the default names in [nextflow.config](nextflow.config) - if file names in the configuration files are altered, the file names should be updated to match.

If you already have any of these data files or reference genomes present locally, e.g. from a prior deployment or other local work, the script will not re-download them, so you can run the script multiple times without worrying about overwriting existing files. This behaviour requires the file names to match exactly (e.g. the reference genome is expected as `ref.fa`).

> **n.b.** Talos has only been tested on GRCh38, and that reference genome is used in a number of places. The download script fetches a gzip-compressed copy of HG38 as part of the setup process, and Talos is not technically supported on other reference genomes, so your variant data will have to be aligned to GRCh38 for the time being. To use the reference genome downloaded by this scripted process, you will need to decompress with `gunzip` and index with `samtools faidx` prior to its use in the workflow.
>
> **n.b.2.** The reference genome is expected to be named `ref.fa` after decompression, and the index file `ref.fa.fai`. If you have an appropriate reference genome already present, you can copy it to the large_files folder and rename it to `ref.fa` (with a matching index) to prevent the script from re-downloading a duplicate reference genome.

## Files

The `gather_files.sh` script downloads the following:

1. A compressed GRCh38 reference genome, in gzipped FASTA format, e.g. from `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/` (see note above about decompressing and indexing)
2. An Echtvar-formatted gnomAD 4.1 reference from `https://zenodo.org/records/22297079`
3. An Ensembl GFF3 file, e.g. `Homo_sapiens.GRCh38.116.gff3.gz` from the [Ensembl FTP site](https://ftp.ensembl.org/pub/release-116/gff3/homo_sapiens).
4. A MANE Summary text file, e.g. `MANE.GRCh38.v1.5.summary.txt.gz` from the [RefSeq MANE FTP site](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5).
5. Predictions for AlphaMissense, e.g. `AlphaMissense_hg38.tsv.gz` from https://zenodo.org/records/8208688.
6. `genes_to_phenotypes.txt` and `hp.obo` from `https://github.com/obophenotype/human-phenotype-ontology/releases`
7. `phenio.db` from https://data.monarchinitiative.org/monarch-kg/latest/phenio.db.gz, then decompressed
8. Raw ClinVar Variant and Submission summaries
