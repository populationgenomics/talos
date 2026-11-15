# **Talos**

[![Docs](https://img.shields.io/badge/docs-populationgenomics.github.io%2Ftalos-blue)](https://populationgenomics.github.io/talos/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
![Test](https://github.com/populationgenomics/automated-interpretation-pipeline/actions/workflows/test.yaml/badge.svg)

> 📖 **Full documentation:** <https://populationgenomics.github.io/talos/>

## **Overview**

As of version 13.0.0, **Talos** has been rebuilt, removing Spark/Hail dependence in favour of a map-reduce framework running directly on sharded VCF representations. This should be an improvement for all users not operating Dataproc/Spark/Hail Batch clusters. The functionality of the tool is unchanged, but hardware requirements will be much reduced.

If you have previously used Talos, see the [Migrating from Talos 1](#migrating-from-talos-1) section below to understand the differences in inputs and outputs for the two versions. Although this repository is the same, `Talos 1` is used here to indicate changes prior to 13.0.0.

---

**Talos** is a scalable, open-source variant prioritisation tool designed to support automated reanalysis of genomic data in rare disease. It identifies **candidate causative variants in known disease genes** by integrating static annotations (e.g. population frequency, predicted consequence) with dynamic knowledge sources such as ClinVar and PanelApp Australia. Talos applies a set of configurable, rule-based logic modules aligned with ACMG/AMP criteria and prioritises variants consistent with expected mode of inheritance and, optionally, patient phenotype.

While Talos can be used for one-off reanalysis of individual families or cohorts, its core design is optimised for **routine, cohort-scale reanalysis**. By comparing current annotations with prior results, Talos highlights **variants that have become reportable due to newly available evidence**—such as new gene–disease or variant–disease relationships—since the last analysis cycle. This enables timely identification of new diagnoses driven by emerging knowledge, while maintaining a low manual review burden.

Talos is specifically intended to identify **variants in established disease genes that are likely to explain the participant’s condition**. It is not designed to detect novel candidate genes or to interpret variants of uncertain significance outside the context of existing clinical knowledge. This focus improves specificity and supports use in diagnostic and research reanalysis workflows.

A full description of the method and its validation in large clinical and research cohorts is available in our publication in Nature Medicine:

[**https://www.nature.com/articles/s41591-026-04477-5**](https://www.nature.com/articles/s41591-026-04477-5)

---

## **When to Use Talos**

Talos is designed to support **automated reanalysis of rare disease cohorts**, enabling identification of **candidate causative variants in known disease genes** based on the latest available evidence. It is best suited for scenarios where:

- You are performing **routine reanalysis** of undiagnosed individuals (e.g. monthly or quarterly)

- You want to detect **variants that have become reportable** due to updates in gene–disease or variant–disease knowledge

- You aim to **minimise the number of variants requiring manual review** optimising for specificity over-sensitivity

- You are working with **exome or genome sequencing data** from previously analysed research or clinical cohorts

- You need a scalable, reproducible pipeline for **family-based or cohort-scale analysis**

> NB. Talos has been explicitly designed for an iterative re-analysis loop, but is effective as a first-pass analysis tool.

Talos is **not currently designed** for:

- Identifying **novel candidate disease genes** or gene discovery

- Analysing **mosaic variants**, or variants outside standard clinical reporting regions

    - Mitochondrial analysis has now been added to Talos, but is limited to ClinVar pathogenic variants only.
    - STR analysis has been added, but currently this is limited to a 'joint-callset' made by combining one or more STRipy JSON results into a VCF format. See [this script](src/talos/scripts/stripy_json_to_vcf.py).

Talos complements existing variant curation workflows by focusing on high-specificity identification of variants that are likely to explain a participant’s condition, based on established gene–disease associations and up-to-date variant-level evidence.

---

## **🚀 Quick Start**

Talos is implemented using **Nextflow**, with all dependencies containerised via Docker. The example workflows can be run locally or on a cluster.

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

Talos requires several large external resources (e.g. reference genome, gnomAD, AlphaMissense, Phenotype data). These are expected in a `large_files` directory. The script [large_files/gather_files.sh](large_files/gather_files.sh) will handle the download of all required resources.

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
- **pedigree**: path to a Pedigree for the cohort, See details [here](docs/Pedigree.md)
- **config**: default available, path to the Talos config - see [example_config.toml](src/talos/example_config.toml) for an example, and the [Configuration README](docs/Configuration.md) for a full breakdown of all config parameters
- **history**: optional, path to previous results
- **ext_ids**: optional, path to ID mapping to present alternate IDs in the HTML report
- **seqr_map**: optional, path to ID mapping to generate hyperlinks to Seqr in the HTML report
- **mito**: optional, path to mitochondrial variants joint-called VCF
- **sv**: optional, path to joint-called SV VCF
- **sv_annotated**: optional, path to an already annotated SV VCF, which skips SV annotation. Mutually exclusive with `sv`, see [SV Annotation](docs/SvAnnotation.md#pre-annotated-sv-vcfs)

The TSV file can contain any number of rows, each representing a distinct Cohort. A parallel Annotation & Talos run will be triggered for each input row, writing to a distinct output folder. An example TSV file has been provided to demonstrate.

#### Input Types 📂

The input TSV uses two columns to locate variant input; `path` and `type`. `path` is the location of the input file or directory. `type` is one of 3 values, **vcf, shards, ss_vcf_dir**

1. **vcf** a single multisample VCF. This will be split into shards and processed in parallel.
2. **shards** a directory of pre-sharded multisample VCF fragments, each shard containing all samples. The annotation and analysis workflows will be parallelised across these shards.
3. **ss_vcf_dir** single-sample VCFs, to be merged in the workflow, then sharded. These are detected using a glob, with the file extension controlled by `params.input_vcf_extension` (defaults to "vcf.bgz")

Outputs are written per cohort: annotation products (annotated VCF shards, a shard manifest, and the annotated SV VCF where relevant) to `{outdir}/{cohort}_annotated`, and analysis results to a dated `{outdir}/{cohort}_analysis_YYYYMMDD`. A proposed input TSV for the next reanalysis cycle, `talos_input_YYYYMMDD.tsv`, is written to the root of `{outdir}`. The output directory should be outside this repository.

The [main.nf](main.nf) workflow runs annotation and analysis together. Annotation is skipped for any cohort whose shard manifest already exists in `{cohort}_annotated`, so the same command is used for the first run and for every reanalysis cycle:

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

## **🔬 Input Validation**
The first step of the Talos Analysis workflow is a module called *StartupChecks*, which runs a number of input validations:

1. Checks a config file is present, and checks all required entries are present and have the correct type
2. Opens an Annotated VCF shard, checking the schema and data types by row sampling and header parsing
3. Parses the Pedigree file and checks that it's well formatted and affected participants are present

This module will either run and complete, or run and fail, printing a collection of all encountered errors. If it fails, you will need to fix the errors before restarting the workflow.

---

## **⚙️ Configuration**

Talos as a Python application is configured through a `TOML` file. This contains all thresholds and parameters for the steps of the Talos workflow. See [`example_config.toml`](src/talos/example_config.toml) as a baseline example, and [Configuration.md](docs/Configuration.md) for extended details on the role of each parameter, and its default value.

Talos as a Nextflow workflow is configured through [`nextflow.config`](nextflow.config) at the repository root, shared by both entrypoints (`preparation.nf`, `main.nf`). It defines paths to annotation resources, container images, and per-process runtime resources; cohort names and per-cohort inputs are supplied through the input TSV instead. [NextflowConfiguration.md](docs/NextflowConfiguration.md) contains a full description of the default values and role in the analysis.

## **📄 Outputs**

Talos produces structured outputs to support both manual review and downstream integration.

- *.json file listing all candidate variants for each proband

  - Includes variant-level and gene-level evidence, inheritance checks, and phenotype match tags

- **HTML reports** summarising results for analysts or clinicians

### **Reanalysis Metadata:**

- first_seen: when the variant was first returned

- evidence_last_updated: when its evidence last changed

Only variants passing configured thresholds and logic modules are returned.

---

## **🔁 Reanalysis Mode**

Talos is designed to support **automated, iterative reanalysis** of undiagnosed cohorts. To do this it reads the results of previous analyses, and integrates them into the latest report. This is currently done by reading in prior analysis results, and incorporating the previous observations with each run. To use this behaviour, use the `history` column in the input TSV - see [History](docs/Reanalysis.md) for more information.

### **How it works:**

1. Run full annotation + prioritisation once

2. In future cycles, keep ClinVar / PanelApp up to date **using the prep workflow**

3. Rerun prioritisation, with the `history` column of the input TSV pointing to a previous run's results, to return all passing variants, with an accurate first-discovered date to enable quick variant triage.

By referencing previous analysis results in each new run, each variant in the output includes:

- first_seen: original detection date

- evidence_last_updated: last evidence update (ClinVar, PanelApp)

> Talos maintains low review burden by allowing users to filter to only variants with newly actionable evidence in each analysis.

---

## **🧬 Phenotype Matching**

Talos supports phenotype-driven filtering using **HPO terms**. See [Pedigree.md](docs/Pedigree.md) for details on how to provide phenotype data in the pedigree file.

### **Matching Strategies:**

- **Patient-to-Gene**: semantic similarity between HPO terms and gene annotations

- **Patient-to-Panel**: PanelApp panels assigned if patient terms match panel HPO tags

- **Cohort-to-Panel**: manually assign panels to all individuals in config

When provided, phenotype terms are used to:

* Build a more accurate gene panel for each analysis, working with PanelApp to match disease-focused panels to HPO terms
* Prioritise variants in the HTML by highlighting variants in genes which are on disease-specific panels, or where the participant and Gene HPO termsets share phenotypic similarities

---

## **🧠 Variant Logic Modules**


Talos prioritises variants using rule-based **logic modules**, each aligned with specific ACMG/AMP evidence criteria.


### **Module Types**

- **Primary**: sufficient to trigger reporting on their own (e.g. ClinVar_PLP)

- **Supporting**: used only as second hits in recessive genes (e.g. AlphaMissense)

### **Standard Modules**

| **Module**          | **Description**                                                                                                                      |
|---------------------|--------------------------------------------------------------------------------------------------------------------------------------|
| ClinVar P/LP        | Pathogenic or Likely Pathogenic by ClinvArbitration                                                                                  |
| ClinVar Recent Gene | P/LP in a PanelApp “new” gene (became Green within the recency window configured by `GeneratePanelData.within_x_months`, default 24) |
| High Impact         | Predicted high-impact protein consequences                                                                                           |
| De Novo             | Confirmed de novo in affected individual                                                                                             |
| PM5                 | Missense in codon with known pathogenic variant                                                                                      |
| LofSV               | Predicted loss-of-function structural variant                                                                                        |
| ClinVar 0-star      | P/LP with 0 gold stars in ClinVar [Supporting category]                                                                              |
| AlphaMissense       | AlphaMissense-predicted pathogenic missense variant  [Supporting category]                                                           |

Each module can be configured through the `.toml` config file (see [Configuration.md](docs/Configuration.md)), up- or down-rating the significance of any individual module depending on the preferences of the deployment host.

---


## **Migrating from Talos 1**

A couple of key points up front:

1. The Pydantic data models being used are either unchanged or have been updated in sync, so output files generated by Talos 1 & 2 should be interchangeable (e.g. Talos1 results as the 'history' file will work)
2. The annotation step needs to be completed in full using the new workflow - the previous stopping point was annotated and reformatted MatrixTables, which are not used from Talos 13.0.0 onwards
3. The Docker image required is completely different and needs to be built fresh, albeit with a much smaller footprint and build time than Talos 1.
4. Several reference files are similar but different - e.g. gnomAD AFs are applied from a new file in the same format, and the Talos 1 data can't be re-used. 

### Migration pathway:

 1. download latest resources
 2. build new docker image
 3. migrate any non-default parameters from the previous `nextflow.config` to this repository's version, e.g. Process/Profile data, or output paths
 4. if you used a custom config TOML file, update any changed keys (`RunHailFiltering` -> `RunSmallFiltering`, `RunHailFilteringSv` -> `RunSvFiltering`). 
 5. run `preparation.nf` workflow
 6. run `main.nf` to complete the annotation through to analysis
 7. on subsequent runs, re-run `main.nf` and annotated shards will be detected and used if available

### Input TSV:

The input TSV file is mostly unchanged from Talos 1, **cohort, path, type, pedigree, and config** are still mandatory columns.

There is one key difference - the 'asset' files (e.g. `nextflow/assets/NO_FILE`) which Talos 1 required as placeholders for missing inputs are no longer required. Instead, any missing columns in the input TSV are automatically interpreted as missing values, and any empty row values are interpreted as missing values for the cohort/row.

e.g. these two config files (as TSVs) are functionally equivalent:

| cohort | path               | type | pedigree              | config               | history  | ext_ids | seqr_map | mito | sv | str |
|--------|--------------------|------|-----------------------|----------------------|----------|---------|----------|------|----|-----|
| CohA   | /path/to/a.vcf.bgz | vcf  | /path/to/pedigree.ped | /path/to/config.toml |          |         |          |      |    |     |

| cohort | path               | type | pedigree              | config               |
|--------|--------------------|------|-----------------------|----------------------|
| CohA   | /path/to/a.vcf.bgz | vcf  | /path/to/pedigree.ped | /path/to/config.toml |

Missing/supplied content can differ per row; this config would provide all mandatory inputs plus the `history` and `mito` file for `CohA`, and only mandatory inputs for `CohB`:

| cohort | path               | type | pedigree              | config               | history         | ext_ids | seqr_map | mito                 | sv | str |
|--------|--------------------|------|-----------------------|----------------------|-----------------|---------|----------|----------------------|----|-----|
| CohA   | /path/to/a.vcf.bgz | vcf  | /path/to/pedigree.ped | /path/to/config.toml | /path/to/a.json |         |          | /path/to/mito.vcf.gz |    |     |
| CohB   | /path/to/b.vcf.bgz | vcf  | /path/to/pedigree.ped | /path/to/config.toml |                 |         |          |                      |    |     |

> **Note**
> the SV reference resources (`svafotate_bed`, `mane_gtf`, `svannotate_noncoding_bed`) must be present for every run, even if no cohort has an `sv` VCF. SV annotation is decided per row: cohorts with an empty or absent `sv` cell skip it.

The optional `str` column takes a bgzipped multi-sample STR VCF, as written by `talos.scripts.stripy_json_to_vcf` from
per-sample STRipy JSON reports. It is not annotated, and is passed straight to `ValidateMOI`. See
[NextflowConfiguration.md](docs/NextflowConfiguration.md#short-tandem-repeats) for details.

### Outputs:

 - `{outdir}/{cohort}_annotated` contains the annotated VCF shards and a manifest file.
 - `{outdir}/{cohort}_analysis_YYYYMMDD` will contain the analysis outputs (HTML/JSON).
 - `{outdir}/talos_input_<runDate>.tsv` is a file containing the same rows as the input for this run, but with the 'history' column pointing to output of the current analysis. This is designed to simplify setting up the next run, so long as the filepaths remain static.

## **📓 Citation**

If you use Talos in your research or clinical workflow, please cite:

> Welland MJ, Ahlquist KD, De Fazio P, et al. _Scalable automated reanalysis of genomic data in research and clinical rare disease cohorts._ Nat Med (2026). https://doi.org/10.1038/s41591-026-04477-5


BibTeX:

```
@article{welland2026talos,
  title     = {Scalable automated reanalysis of genomic data in research and clinical rare disease cohorts},
  author    = {Welland, Matthew J and Ahlquist, KD and De Fazio, Paul and Austin-Tse, Christina and Pais, Lynn and Wedd, Laura and Bryen, Samantha and Rius, Rocio and Franklin, Michael and Hall, Giles and et al.},
  journal   = {Nature Medicine},
  year      = {2026},
  doi       = {10.1038/s41591-026-04477-5},
  url       = {https://doi.org/10.1038/s41591-026-04477-5},
}
```
