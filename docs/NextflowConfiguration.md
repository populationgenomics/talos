# Nextflow Configuration

This README documents the parameters used in the nextflow configuration file. Some default values for configuration variables are initially set up to facilitate the demonstration test-run.

| **Parameter**           | **Description**                                                                                                                                                                     |
|-------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `large_files`           | Directory containing large external resources (e.g. gnomAD, MANE, AlphaMissense). See [large_files README](../large_files/README.md) for instructions on how to source this content |
| `processed_annotations` | Output directory for all processed/reformatted data. Cohort/callset independent, to be used across multiple runs.                                                                   |
| `outputDir`             | Defaults to `processed_annotations`, required for the 'workflow outputs' NextFlow functionality. Override for Annotation/Talos workflow using `-output-dir <path>`                  |
| `input_vcf_extension`   | Linked to `shards`. If a sharded VCF is provided, this argument is the file extension to glob for. By default this is `vcf.bgz`, but could be altered to `vcf.gz`.                  |
| `alphamissense_tsv`     | Path in `large_files` to the AlphaMissense raw data, reformatted to be used in Hail annotation                                                                                      |
| `alphamissense_zip`     | Path in `processed_annotations` to the reformatted AlphaMissense data, ready to be used by Echtvar for annotation                                                                   |
| `ensembl_bed`           | Path in `processed_annotations` to the Ensembl BED file, one line per gene                                                                                                          |
| `ensembl_merged_bed`    | Path in `processed_annotations` to the Ensembl BED file, overlapping regions merged                                                                                                 |
| `mane_json`             | Path in `processed_annotations` to pre-processed MANE data                                                                                                                          |
| `ensembl_gff`           | Path in `large_files` to the Ensembl GFF3 file                                                                                                                                      |
| `gnomad_zip`            | Path in `large_files` to the echtvar annotation file from [Zenodo](https://zenodo.org/records/15222100)                                                                             |
| `mane`                  | Path in `large_files` to the MANE transcript data                                                                                                                                   |
| `ref_genome`            | Path in `large_files` to the reference genome FASTA file                                                                                                                            |
| `hpo`                   | Path in `large_files` to the HPO ontology file in `.obo` form                                                                                                                       |
| `gen2phen`              | Path in `large_files` to the Geno~Phenotype databset in `.txt` form                                                                                                                 |
| `phenio_db`             | Path in `large_files` to the Monarch Phenotype DB in `.db.gz` form                                                                                                                  |
| `submission_summary`    | FTP URL for weekly ClinVar submission dump                                                                                                                                          |
| `variant_summary`       | FTP URL for weekly ClinVar variant dump                                                                                                                                             |
| `clinvar_blacklist`     | A string, containing quoted, space-delimited entries for all ClinVar submitter sites to blacklist (ignore).                                                                         |
| `container`             | Docker image to use                                                                                                                                                                 |
| `docker.enabled`        | Parameter for Nextflow to enable usage of Docker                                                                                                                                    |
