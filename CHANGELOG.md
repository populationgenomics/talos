# Changelog

All notable changes to this project will be documented in this file.

Suggested headings per release (as appropriate) are:

* `Added` for new features.
* `Changed` for changes in existing functionality.
* `Deprecated` for soon-to-be removed features.
* `Removed` for now removed features.
* `Fixed` for any bug fixes.
* `Security` in case of vulnerabilities.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

<!--changelog-start-->
<!--latest-start-->

[13.0.0]

## Added

  * Optional `sv_annotated` input TSV column: an already annotated SV VCF (e.g. GATK-SV output) which bypasses SV annotation and goes straight to SV filtering. Setting both `sv` and `sv_annotated` on one row is an error
  * Small-variant and `sv` VCFs are subset to the pedigree's samples during annotation, dropping variants no pedigree sample carries, so a dataset-wide callset is annotated and stored at cohort size. VCFs already matching the pedigree are unchanged. Pedigree samples missing from a VCF are a warning, no shared samples an error. See `docs/Pedigree.md`
  * `params.recalculate_ac_af_an` (default `true`). Set `false` to keep the INFO/AC, AF and AN already in the input VCF, e.g. from a larger callset, instead of recalculating them from the cohort's genotypes. Pedigree subsetting then leaves them alone too

## CPG only

  * `talos.cpg_internal.input_tsv_from_metamist` builds input TSV rows for Metamist cohorts, generating each cohort's pedigree, external ID map, seqr map and merged config. It needs the new `cpg` extra (`metamist`)
  * `cpg_main.nf`, a second entrypoint over `TALOS_PIPELINE`. It takes Metamist cohort IDs (`--cohort_ids`, `--access_level`, `--sequencing_type`, `--talos_config`) instead of an input TSV, builds the cohort rows with `InputTsvFromMetamist`, and publishes each cohort to `gs://cpg-<dataset>-{main|test}/<sequencing_type>/talos/<cohort>/{annotated,analysis_YYYYMMDD}/`. See `docs/CpgDeployment.md`
  * `cpg_main.nf --metamist` registers each cohort's published results JSON (`aip-results`, the next run's history) and HTML report (`aip-report`) in Metamist, against the Metamist cohort. Off by default

## Changed

  * Outputs are published with process `publishDir` instead of workflow outputs (the `output {}` block), so each cohort can be given its own output location. The entrypoint decides each cohort's annotated and analysis directories and passes them down to the processes that publish. The output layout is unchanged
  * `main.nf` takes its output root from the required `--outdir` param. `-output-dir` / `outputDir` is no longer read
  * `preparation.nf` publishes each newly generated resource directly to `params.processed_annotations`, and no longer needs an output directory
  * The annotation → SV annotation → analysis wiring moved out of `main.nf` into a shared `TALOS_PIPELINE` workflow (`nextflow/talos_pipeline.nf`), which takes one row per cohort with its output directories. `main.nf` now only validates params, builds those rows from the input TSV, and writes the next-cycle TSV. Process names gain a `TALOS_PIPELINE:` prefix, so `-resume` won't reuse tasks cached by earlier runs
  * SV annotation is decided per input TSV row instead of by whether the TSV has an `sv` column. The SV reference data (`svafotate_bed`, `mane_gtf`, `svannotate_noncoding_bed`) is now required for every run
  * Some TOML config sections renamed from their Hail versions (e.g. RunHailFiltering -> RunSmallFiltering)

## Fixed

  * Single-sample VCF inputs (`ss_vcf_dir`) are now split and left-aligned per sample before merging, instead of after (#729). When an allele was part of a multi-allelic record in one sample but a separate record in another, the merge split it across two rows, each with only some samples' genotypes. That produced false de novo calls and dropped alleles. Cohorts annotated from `ss_vcf_dir` before this fix keep the bad rows until their `${cohort}_annotated/` directory is deleted and they are re-annotated. A variant written twice within one sample's VCF is now also deduplicated
  * The default `errorStrategy` in `nextflow.config` read `process.maxRetries`, which is undefined when the closure runs. Every task failure threw instead of retrying. It now reads `task.maxRetries`, so failed tasks are retried up to 3 times with backoff

## Removed

  * `params.annotated_dir`. Annotation products are always read from, and written to, `<outdir>/<cohort>_annotated`
  * The startup check comparing the output directory with `params.processed_annotations`
  * `workflow.output.mode` and `outputDir` from `nextflow.config`

[0.6.0] - 2026-10

## Added

  * PanelApp CNV regions tagged `regulatory region` are downloaded and carried through to the parsed PanelApp data
  * `run_sv_filtering.py` labels rare DEL/DUP CNVs matching a region's type (`cnv_loss`/`cnv_gain`/`cnv_both`) overlapping at least the region's required % of the region (CNV size doesn't matter) as `categorybooleanregulatorycnv`. Each match gets its own row, using the region name as the gene ID, so region matches are retained without any gene annotation, and are never paired with gene variants as comp-hets
  * The report shows region name, location, type, observed overlap and MOI for region matches
  * CNVs fully contained in a region must still cover at least `RunSvFiltering.regulatory_region_percentage` (default 10) % of the region to match
  * `generate_sv_test_data.py` adds a positive (DUP) and negative (DEL) control for the IHH upstream regulatory region, on panel 93, now forced in `nextflow/inputs/config.toml`

## Changed

  * SV input no longer needs every GATK-SV field. `AC`, `AF`, `AN`, `N_HET`, `N_HOMALT`, `MALE_AF`/`FEMALE_AF` and `SVLEN` are optional; a missing frequency field logs a warning and its filter is skipped, except the per-sex AF filter, which falls back to overall `AF`
  * `prepare_sv_for_svannotate.py` backfills a missing `SVLEN` as `END - POS` where `END > POS`

[0.5.0] - 2026-09

## Added

  * STR data can be supplied per cohort through a new optional `str` column in the input TSV, and is passed to ValidateMOI
  * `nextflow/inputs/generate_str_test_data.py` simulates an STR VCF for the test trio, used by `test_sv.tsv`

## Changed

  * The way STRs are handled has been changed. PanelApp now contains STRs curated into each phenotype-specific panel
  * STRs are detected on each individual panel, instead of all being sourced from the 'repeat disorders' holding panel
  * By default, STR analysis is only done for phenotype-matched individuals

## Fixed

  * GATK SVAnnotate only works with symbolic alleles. `SortCpxIntervals` is renamed `PrepareSvForSvAnnotate`, and now also removes any SV record with a non-symbolic ALT before annotation, logging each removed record. Update any `withName: SortCpxIntervals` resource blocks accordingly

[0.4.0] - 2026-09

### Added

* Panel confidence (Red/Amber/Green) is now recorded in the history, and if a variant is newly promoted to Green, this is identified in the report
* The Report has Emojis next to the date if variant evidence changes:
    ➕ - Seen previously, but has additional categories in this round (briefly 🆕, now a new emoji which better represents 'addition' instead of a totally new variant)
    💫 - ClinVar star count has increased
    🎯 - Newly phenotype-matched
    🚦 - Newly Green in PanelApp

[0.3.1] - 2026-09

### Changed

* The MakeScatterRegions module only detects and scatters across chr1-22,X,Y, ignoring chrM and decoy contigs. There are no PanelApp genes on decoy contigs, and chrM has a separate annotation flow

[0.3.0]

### Changed

* The annotation workflow now runs echtvar and `bcftools csq` on a sites-only copy of each normalised shard (all samples dropped), then lifts the resulting INFO fields back onto the full-width shard with `bcftools annotate`. Only `bcftools norm`/`+fill-tags` and the final lift touch the genotype columns, so the per-shard cost no longer scales with cohort width for the annotation tools themselves. `AnnotateWithEchtvar` and `AnnotateCsqWithBcftools` are replaced by a single `AnnotateShard` process - update any `withName:` resource blocks accordingly.
* The gnomAD AF < 0.05 pre-filter has moved from echtvar's `-i` expression to a `bcftools view -e 'INFO/gnomad_AF >= 0.05'` on the annotated full-width shard. Variants absent from gnomAD (`gnomad_AF=0`) are retained exactly as before.
* Annotated shards are now published with a `.tbi` index alongside. Existing `${cohort}_annotated` directories (without indexes) remain valid - the manifest still names the shards to reuse.
* The gnomAD 4.1 annotation source has been changed - the new Zenodo record is https://zenodo.org/records/22297079, and it is referenced by a new name in the config file to differentiate the two versions. The new version of the data contains fewer annotations, but spans every base in the gnomAD dataset, which prevents issues where frequency annotations could have been missing due to over-zealous region filtering when generating the original record.

* `NormaliseVcf` writes its full-width intermediate at the lightest BCF compression level with threaded deflate, and `AnnotateShard` compresses its output across `task.cpus` (now 2). Deflate was ~70% of the per-shard runtime.

### Removed

* Annotated shards no longer carry a per-sample `FORMAT/BCSQ` bitmask. Nothing in Talos read it; `INFO/BCSQ` is unchanged.
* `FORMAT/PL` is dropped from every shard during normalisation, and so is absent from the annotated and labelled VCFs. Talos never read it (only GT, GQ, DP, AD and the phasing fields), and it was the widest per-sample field - published shards shrink by roughly a fifth.

[0.2.0] - 2026-08

### Added

* Adds a new category for BND variants occurring in Exons

[0.1.0] - 2026-08

* Every `main.nf` run now writes `talos_input_YYYYMMDD.tsv` to the root of the output directory - a copy of the run's input TSV with each cohort's `history` cell repointed at the results JSON that run published, ready to be reviewed and used as the `--input_tsv` for the next reanalysis cycle.
* This aims to bridge the gap in the current implementation, where inputs and analysis outputs need to be manually updated for future cycles. Feel free to use this or ignore it :)

[0.0.2] - 2026-08

### Changed

* Really reduces the per job resource requirements. This may be insufficient in HUGE cohorts, but has been successful in mid-sized (50-genome) testing. Users are invited to customise based on local requirements.

### Fixed

* Fixes the mechanism NextFlow uses to find adjacent-files, and carry out file globbing. This was fine for local deployments, but buggy in cloud environments due to String trimming.

[0.0.1] - 2026-08

### Added

* Welcome to the new repository! This is Talos 2.0, designed to be easier to deploy outside the CPG and Hail Batch-using sites. Instead of using Hail/Spark as a processing framework, this reimagined version uses a more basic map/reduce framework, scattering the input data as sharded VCFs, annotating and filtering in pieces, and re-forming a final VCF input from the minimal set of relevant variants. No Hail runtime, no MatrixTables, no Spark cluster. Just VCFs. This feels like a good compromise given the considerations - this needs to remain site-agnostic, making few assumptions about available infrastructure, but needs to run relatively quickly across cohorts of all sizes. This will still present scaling issues at super-high cohort sizes, but it represents dramatic improvements over the original Talos implementation at non-Spark enabled sites.

### Changes from the original Talos

* The labelling/filtering phase which runs on the annotated data has been rewritten as a loop over an annotated VCF fragment, instead of requiring data to be read into a MatrixTable and processed in a single highly-resourced runtime.
* The annotation workflow no longer physically splits the input VCF (`SplitVcf` is gone - it was a serial pass over the whole callset). Scatter regions are now derived from the VCF index (`MakeScatterRegions`), with each `NormaliseVcf` task reading only its region from the input; `params.vcf_split_n` is now the target records per region, and `0` still disables the scatter.
* `MakeScatterRegions` now stages only the `.tbi` - `bcftools index -s` reads contig, length and record count from the index alone, so the whole-callset VCF is no longer localised into that task.
* Outputs are split per cohort: annotation products (annotated shards, shard manifest, annotated SV VCF) publish to `{cohort}_annotated/`, and analysis results to a dated `{cohort}_analysis_YYYYMMDD/`. The `{cohort}_annotated` directory is reusable across reanalysis cycles and safe to delete to force re-annotation. **Migration**: existing output directories keep annotation products in `{cohort}_outputs/` with no manifest - re-run `main.nf` once to regenerate into the new layout.
* Annotation writes a per-cohort `{cohort}_manifest.json` naming exactly the shards it produced, and is now idempotent: cohorts whose manifest already exists are reused (shards read from the manifest, never a directory glob), so stale shards from an earlier annotation run (e.g. a different `vcf_split_n`) can never be double-counted, and a reanalysis cycle is simply a re-run of `main.nf`. A new optional `params.annotated_dir` lets a run read annotation products from a different location than it publishes results to.
* VCF ingestion in ValidateMOI no longer iterates the full cohort three times per variant - depths, alt depths and AB ratios are sliced from the genotype arrays at carrier positions only, and phase data is read for carriers only. Per-variant config lookups (exomiser rank, de novo thresholds, CSQ field names) are now resolved once per process.

### Removed

* `talos_only.nf` - redundant now that annotation is idempotent; re-run `main.nf` for reanalysis cycles instead.

<!--latest-end-->
<!--changelog-end-->
