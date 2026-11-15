# Getting Started

This guide walks you through preparing the environment, downloading the required reference data, and running your first Talos analysis.

Talos is implemented using **Nextflow**, with all dependencies containerised via Docker. The example workflows can be run locally or on a cluster.

## Resource Configuration

This re-worked version of Talos removes the need to operate a large Spark cluster, instead running the majority of filter and labelling operations
in streamed processes, running on fragments of the whole dataset, using efficient third-party tools. Instead of requiring large VMs, this changes
to a map-reduce strategy to solve large data issues. As such, the need for large VMs to run individual steps is greatly reduced. You may still
benefit from increasing memory, but it shouldn't be the hard requirement it was previously.

The standard process resourcing in the `nextflow.config` file should be sufficient for most users, with only a few steps
given extra resource allocations:

| Stage                   | Workflow      | Suggestion           | Reasoning                                                                                                                    |
|:------------------------|:--------------|:---------------------|:-----------------------------------------------------------------------------------------------------------------------------|
| AnnotateSvWithSvafotate | SV Annotation | 16.GB memory, 4 CPUs | SV Annotation may not be utilised by most users, but when required it involves loading a large reference object into memory. |
| ValidateMOI             | Talos         | High memory          | Runs across the whole dataset, forming DataClass instances from each variant. In large cohorts this can run hot.             |

Honorable mentions:

`params.vcf_split_n` (default: `2500000`) is the target number of records per scatter region. Regions are derived from the VCF index (no physical split pass); each annotation task then reads only its region from the input VCF. For datasets with very high sample counts, reducing this value will create more regions, better for the parallelised annotation and filtering workflow. Setting it to `0` disables the scatter entirely.

---

## Workflow

There are two primary workflows:

* `preparation.nf`: downloads and formats data in preparation for Talos runs.
* `main.nf`: imports and executes the two sub-workflows which comprise the Talos runs - Annotation, and Analysis. Annotated outputs are idempotent, and can be re-used in future runs.

There are 3 main file paths you need to provide. The default locations for these are in the project root, but for real analyses you will want them to live outside the project folder:

* `params.large_files` or `--large_files`: the directory containing large raw annotation sources. This defaults to `large_files` in the project root for the test workflow to run, but for regular use this should be in a storage drive accessible to your HPC/Cloud/local compute environment. All workflows need to reference these files.
* `params.processed_annotations` or `--processed_annotations`: the outputs of the preparation sub-workflow, which processes raw 'large files' into annotation-ready resources
* `params.outdir` or `--outdir`: the root folder/bucket for the main workflow's per-cohort outputs. Required by `main.nf`, unused by `preparation.nf`, which always writes to `params.processed_annotations`.

### **1. Install Requirements**

- [Nextflow](https://www.nextflow.io/docs/latest/install.html)

- Docker

To build the Docker image:

```
docker build -t talos:13.0.0 .
```

### **2. Download Annotation Resources**

Talos requires several large external resources (e.g. reference genome, gnomAD, AlphaMissense, Phenotype data). These are expected in a `large_files` directory. The script [large_files/gather_files.sh](https://github.com/populationgenomics/talos/blob/main/large_files/gather_files.sh) will handle the download of all required resources.

### **3. Run Preparation Workflow**

The preparations workflow transforms raw data into usable annotation sources. Most of these only need to be prepared once (MANE, GFF3, AlphaMissense), but the ClinVar and PanelApp data need to be downloaded once per month to keep the workflow up to date.

```bash
nextflow \
    -c nextflow.config \
    run preparation.nf \
    [--processed_annotations <processed_annotations_path>] \
    [--large_files <large_files_path>]
```

### **4. Run Annotation & Talos Combined Workflow**

> Cohort-specific inputs for the Talos workflow are now provided in a single file, `--input_tsv`, instead of using several separate parameters.

The inputs for the Talos workflow are:
- **cohort**: a collective name to identify the input/results, used in output directory and file naming
- **path**: path to the Cohort's input data (VCF)
- **type**: type of the input data, see below
- **pedigree**: path to a Pedigree for the cohort, See details [here](Pedigree.md)
- **config**: default available, path to the Talos config - see [example_config.toml](https://github.com/populationgenomics/talos/blob/main/src/talos/example_config.toml) for an example, and the [Configuration README](Configuration.md) for a full breakdown of all config parameters
- **history**: optional, path to previous results
- **ext_ids**: optional, path to ID mapping to present alternate IDs in the HTML report
- **seqr_map**: optional, path to ID mapping to generate hyperlinks to Seqr in the HTML report
- **mito**: optional, path to mitochondrial variants joint-called VCF
- **sv**: optional, path to joint-called SV VCF
- **sv_annotated**: optional, path to an already annotated SV VCF, which skips SV annotation. Mutually exclusive with `sv`, see [SV Annotation](SvAnnotation.md#pre-annotated-sv-vcfs)

The TSV file can contain any number of rows, each representing a distinct Cohort. A parallel Annotation & Talos run will be triggered for each input row, writing to a distinct output folder. An example TSV file has been provided to demonstrate.

#### Input Types 📂

The input TSV uses two columns to locate variant input; `path` and `type`. `path` is the location of the input file or directory. `type` is one of 3 values, **vcf, shards, ss_vcf_dir**

1. **vcf** a single multisample VCF. This will be split into shards and processed in parallel.
2. **shards** a directory of pre-sharded multisample VCF fragments, each shard containing all samples. The annotation and analysis workflows will be parallelised across these shards.
3. **ss_vcf_dir** single-sample VCFs, to be merged in the workflow, then sharded. These are detected using a glob, with the file extension controlled by `params.input_vcf_extension` (defaults to "vcf.bgz")

Outputs are written per cohort: annotation products (annotated VCF shards, a shard manifest, and the annotated SV VCF where relevant) to `{outdir}/{cohort}_annotated`, and analysis results to a dated `{outdir}/{cohort}_analysis_YYYYMMDD`. A proposed input TSV for the next reanalysis cycle, `talos_input_YYYYMMDD.tsv`, is written to the root of `{outdir}`. The output directory should be outside this repository.

The [main.nf](https://github.com/populationgenomics/talos/blob/main/main.nf) workflow runs annotation and analysis together. Annotation is skipped for any cohort whose shard manifest already exists in `{cohort}_annotated`, so the same command is used for the first run and for every reanalysis cycle:

```bash
nextflow \
  -c nextflow.config \
  run main.nf \
  --input_tsv nextflow/inputs/test.tsv \
  --processed_annotations <processed_annotations_path> \
  --large_files <large_files_path> \
  --outdir <path_to_output_dir>
```

>**For best results we advise repeating the Talos workflow on a regular cadence** — just re-run the command above; delete a cohort's `{cohort}_annotated` directory to force re-annotation (e.g. when MANE data is updated).

---

Outputs are split into two directories per cohort:

* `{outdir}/{cohort}_annotated` — annotation products: the annotated VCF shards, a
  `{cohort}_manifest.json` naming exactly the shards produced, and the annotated SV VCF where SV data was
  supplied. The annotation sub-workflow only needs to be run once per dataset — these files are reused for every
  subsequent reanalysis cycle, and the whole directory is safe to delete to force clean re-annotation without
  touching any analysis results.
* `{outdir}/{cohort}_analysis_YYYYMMDD` — analysis results (report HTML, results JSON, PanelApp
  data, labelled VCFs), one directory per analysis date.

For subsequent cycles, run the same `main.nf` command again: annotation is skipped for every cohort whose
shard manifest is already present, and only the analysis stages re-run.

!!! tip "Reanalysis cadence"
    For best results, re-run the Talos workflow on a regular cadence — monthly or quarterly is typical. See the [Reanalysis](Reanalysis.md) page for how prior results are folded into new runs.

---

## 6. Check the outputs

A successful run produces:

- A `*.json` file listing all candidate variants for each proband, with variant- and gene-level evidence, inheritance checks, and phenotype tags.
- An optional **HTML report** for analyst/clinician review.
- An optional **simplified TSV** for Seqr ingestion via `MinimiseOutputForSeqr`.

Each variant carries reanalysis metadata (`first_seen`, `evidence_last_updated`) so downstream consumers can filter to newly actionable findings.

---

## What happens during startup

The first step of every Talos run is a `StartupChecks` module that validates inputs before any analysis begins:

1. Confirms a config file is present and that all required entries exist with the correct types.
2. Opens the VCF and checks the schema and data types, checking for overlap between the VCF samples and pedigree.
3. Parses the pedigree file, validating its format and that affected participants are present.
4. Checks the ClinVar data, ensuring it is recent and has sufficient entries.

If startup fails, Talos prints a collected list of all encountered errors. Fix these before restarting the workflow.

---

## Next steps

- Tune the rule-based logic via [Talos Configuration](Configuration.md).
- Tune Nextflow runtime and resource settings via [Nextflow Configuration](NextflowConfiguration.md).
- Read [Features](features.md) for an overview of the logic modules and reanalysis behaviour.
