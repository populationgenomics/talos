# Nextflow Configuration

This README documents the parameters used in the nextflow configuration file. Some default values for configuration variables are initially set up to facilitate the demonstration test-run.

| **Parameter**           | **Description**                                                                                                                                                                                                                        |
|-------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `large_files`           | Directory containing large external resources (e.g. gnomAD, MANE, AlphaMissense). See [large_files README](https://github.com/populationgenomics/talos/blob/main/large_files/README.md) for instructions on how to source this content |
| `processed_annotations` | Output directory for all processed/reformatted data. Cohort/callset independent, to be used across multiple runs.                                                                                                                      |
| `outdir`                | Required by `main.nf`. Root for per-cohort outputs: `<outdir>/<cohort>_annotated` (annotation products, reused by later runs — annotation is skipped for cohorts whose manifest exists there) and `<outdir>/<cohort>_analysis_YYYYMMDD` (results). |
| `input_vcf_extension`   | Linked to `shards`. If a sharded VCF is provided, this argument is the file extension to glob for. By default this is `vcf.bgz`, but could be altered to `vcf.gz`.                                                                     |
| `recalculate_ac_af_an`  | Default `true`: INFO/AC, AF and AN are recalculated from the cohort's genotypes (after subsetting to the pedigree) with `bcftools +fill-tags`. Set `false` to keep the values already in the input VCF, e.g. counts from a larger callset; AC and AN must then be present. Applied during annotation, so existing `<outdir>/<cohort>_annotated` shards keep whichever setting they were made with. |
| `alphamissense_tsv`     | Path in `large_files` to the AlphaMissense raw data, reformatted to be used in Hail annotation                                                                                                                                         |
| `alphamissense_zip`     | Path in `processed_annotations` to the reformatted AlphaMissense data, ready to be used by Echtvar for annotation                                                                                                                      |
| `ensembl_bed`           | Path in `processed_annotations` to the Ensembl BED file, one line per gene                                                                                                                                                             |
| `ensembl_merged_bed`    | Path in `processed_annotations` to the Ensembl BED file, overlapping regions merged                                                                                                                                                    |
| `mane_json`             | Path in `processed_annotations` to pre-processed MANE data                                                                                                                                                                             |
| `mitimpact_zip`         | Path in `processed_annotations` to the reformatted MitImpact reference data                                                                                                                                                            |
| `mitotip_zip`           | Path in `processed_annotations` to the reformatted MitoTip reference data                                                                                                                                                              |
| `napogee_zip`           | Path in `processed_annotations` to the reformatted nAPOGEE reference data                                                                                                                                                              |
| `ensembl_gff`           | Path in `large_files` to the Ensembl GFF3 file                                                                                                                                                                                         |
| `gnomad_zip`            | Path in `large_files` to the echtvar annotation file from [Zenodo](https://zenodo.org/records/15222100)                                                                                                                                |
| `mane`                  | Path in `large_files` to the MANE transcript data                                                                                                                                                                                      |
| `ref_genome`            | Path in `large_files` to the reference genome FASTA file                                                                                                                                                                               |
| `hpo`                   | Path in `large_files` to the HPO ontology file in `.obo` form                                                                                                                                                                          |
| `gen2phen`              | Path in `large_files` to the Geno~Phenotype databset in `.txt` form                                                                                                                                                                    |
| `phenio_db`             | Path in `large_files` to the Monarch Phenotype DB in `.db.gz` form                                                                                                                                                                     |
| `mitimpact_tsv`         | Path in `large_files` to the raw Mitimpact file                                                                                                                                                                                        |
| `mitotip_tsv`           | Path in `large_files` to the raw MitoTip file                                                                                                                                                                                          |
| `napogee_tsv`           | Path in `large_files` to the raw Mitochondrial nAPOGEE file                                                                                                                                                                            |
| `submission_summary`    | FTP URL for weekly ClinVar submission dump                                                                                                                                                                                             |
| `variant_summary`       | FTP URL for weekly ClinVar variant dump                                                                                                                                                                                                |
| `clinvar_blacklist`     | A string, containing quoted, space-delimited entries for all ClinVar submitter sites to blacklist (ignore).                                                                                                                            |
| `container`             | Docker image to use                                                                                                                                                                                                                    |
| `docker.enabled`        | Parameter for Nextflow to enable usage of Docker                                                                                                                                                                                       |

## CPG entrypoint

Only `cpg_main.nf` reads these, and there all but `metamist` and `metamist_wait_minutes` are required. `main.nf` ignores them. See
[CPG Deployment](CpgDeployment.md).

| **Parameter**     | **Description**                                                                                                                  |
|-------------------|----------------------------------------------------------------------------------------------------------------------------------|
| `cohort_ids`      | Metamist cohort IDs, as a list or a comma-delimited string                                                                       |
| `access_level`    | `test`, `standard` or `full`. `test` publishes to `gs://cpg-<dataset>-test`, the others to `gs://cpg-<dataset>-main`            |
| `sequencing_type` | `genome` or `exome`                                                                                                              |
| `talos_config`    | The CPG `talos.toml`. Its `[cohorts.<dataset>]` sections are merged into each cohort's config                                    |
| `metamist`        | Optional, default `false`. Register each cohort's results JSON (`aip-results`) and HTML report (`aip-report`) in Metamist      |
| `metamist_wait_minutes` | Optional, default `30`. How long registration waits for each published file                                              |

## Structural Variant annotation

These parameters are used by the SV annotation workflow. Its reference checks run for every run, so the
`svafotate_bed`, `mane_gtf` and `svannotate_noncoding_bed` files must exist even if no cohort has an `sv` VCF.

| **Parameter**              | **Description**                                                                                     |
|----------------------------|-----------------------------------------------------------------------------------------------------|
| `svafotate_bed`            | Path in `large_files` to the SVAFotate population-frequency BED. **Use as downloaded** — see below  |
| `mane_gtf`                 | Path in `large_files` to the MANE GTF, used by GATK SVAnnotate for gene consequences                |
| `svannotate_noncoding_bed` | Path in `large_files` to the GATK-SV non-coding elements BED                                        |
| `ref_dict`                 | Path in `processed_annotations` to the sequence dictionary. Generated on the first SV run if absent |
| `sv_overlap_fraction`      | Reciprocal overlap fraction required to match an SV against the gnomAD reference. Defaults to `0.5` |
| `gatk_container`           | Docker image for GATK SVAnnotate. Pulled from the public registry, not built locally                |

### Providing the input

The joint-called SV VCF is supplied per cohort, as an `sv` column in the input TSV, alongside the existing
`mito` column. A bgzipped VCF with a matching `.tbi` is expected. Cohorts with no SV data simply leave the
column empty.

A full GATK-SV callset carries `AC`, `AF`, `AN`, `N_HET`, `N_HOMALT`, array-typed `MALE_AF`/`FEMALE_AF` (or
`AF_MALE`/`AF_FEMALE`) and `SVLEN`; a bare gCNV or caller-native VCF may not. All are optional. The annotation
chain doesn't write the frequency fields, so when they are absent `RunSvFiltering` logs a warning and skips the
filter that depends on them. Without `MALE_AF`/`FEMALE_AF` the per-sex callset AF filter falls back to overall `AF`
(skipped if that is absent too), and the callset AC/AF filter in `ValidateMOI` is skipped without `AC`/`AF`. Expect more, commoner SVs in the report. A missing `SVLEN` is
backfilled as `END - POS` before SVAnnotate, where `END` is present and greater than `POS`.

### Trying it out

`nextflow/inputs/test_sv.tsv` is the SNV test input plus an `sv` column, pointing at a small adversarial SV
VCF (it also has an `str` column, see [Short Tandem Repeats](#short-tandem-repeats)). `nextflow/inputs/test.tsv` deliberately has no `sv` column, so the default test run never requires the
478 MB of SV reference data. Regenerate the test VCF with:

```bash
python nextflow/inputs/generate_sv_test_data.py
```

### Re-runs are skipped automatically

If `<outdir>/<cohort>_annotated/<cohort>_sv_annotated.vcf.bgz` already exists, that cohort's annotation is skipped
entirely and the existing file is reused. Annotation is expensive, so delete that file to force re-annotation.

### Do not add `chr` prefixes to the SVAFotate BED

The SVAFotate BED uses Ensembl-style contig names (`1`, not `chr1`), which looks inconsistent with the rest of
Talos. It is correct. SVAFotate strips the `chr` prefix from the query VCF but not from this reference file, so
adding prefixes causes **every** variant to be annotated with a population frequency of zero, with no error
raised — making an entire callset appear rare. Use the file exactly as downloaded.

## Short Tandem Repeats

STR data needs no Nextflow parameters or reference files. Supply a multi-sample STR VCF per cohort in the `str` column
of the input TSV. Leave the column empty, or leave it out altogether, for cohorts with no STR data. The VCF is not
annotated or filtered first; it goes straight to `ValidateMOI`, alongside the labelled small-variant, SV, and mito
VCFs.

### Providing the input

The VCF has to be bgzipped. `ValidateMOI` indexes it before reading. Build it from per-sample STRipy JSON reports with:

```bash
python -m talos.scripts.stripy_json_to_vcf \
    --json sample1.json sample2.json ... \
    --output cohort_str.vcf.bgz \
    --mapping GRCh38_symbol_to_ensg.json
```

`--mapping` is the symbol-to-ENSG lookup written to `processed_annotations` by the preparation workflow. It is needed
because STRs are matched to PanelApp by gene ID, and STRipy only reports a locus name.

Talos reads these fields:

| Field                    | Use                                                                                    |
|--------------------------|----------------------------------------------------------------------------------------|
| `INFO/GENE`              | ENSG ID - matched against the PanelApp STR entities on this cohort's panels            |
| `INFO/LOCUS`             | STRipy locus ID, e.g. `ARX_1` - checked against `ValidateMOI.noisy_strs` in the config |
| `FORMAT/GT`              | `1` marks an allele in STRipy's pathogenic range, `0` any other allele                 |
| `FORMAT/REPCN`           | repeat count for each allele, shown in the report                                      |
| `FORMAT/DISEASE_DETAILS` | per-disease normal/intermediate ranges and pathogenic threshold                        |

The only STR calls that can be reported are those at PanelApp STR entities that pass the cohort's `confidence_level`
on one of its panels. By default only phenotype-matched and forced panels count, not the `default_panel`. See
`pheno_match_strs` in [Configuration.md](Configuration.md).

### Trying it out

`nextflow/inputs/test_sv.tsv` includes a `str` column for the `premerged` cohort. It points to a small simulated STR
VCF, whose cases and expected outcomes are listed in the docstring of the generator. The test config sets
`pheno_match_strs = false`, as none of the test proband's matched panels carry an STR. Regenerate the VCF with:

```bash
uv run python nextflow/inputs/generate_str_test_data.py
```
