# **Talos**

[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
![Test](https://github.com/populationgenomics/automated-interpretation-pipeline/actions/workflows/test.yaml/badge.svg)

## **Overview**

**Talos** is a scalable, open-source variant prioritisation tool designed to support automated reanalysis of genomic data in rare disease. It identifies **candidate causative variants in known disease genes** by integrating static annotations (e.g. population frequency, predicted consequence) with dynamic knowledge sources such as ClinVar and PanelApp Australia. Talos applies a set of configurable, rule-based logic modules aligned with ACMG/AMP criteria and prioritises variants consistent with expected mode of inheritance and, optionally, patient phenotype.


While Talos can be used for one-off reanalysis of individual families or cohorts, its core design is optimised for **routine, cohort-scale reanalysis**. By comparing current annotations with prior results, Talos highlights **variants that have become reportable due to newly available evidence**—such as new gene–disease or variant–disease relationships—since the last analysis cycle. This enables timely identification of new diagnoses driven by emerging knowledge, while maintaining a low manual review burden.


Talos is specifically intended to identify **variants in established disease genes that are likely to explain the participant’s condition**. It is not designed to detect novel candidate genes or to interpret variants of uncertain significance outside the context of existing clinical knowledge. This focus improves specificity and supports use in diagnostic and research reanalysis workflows.


A full description of the method and its validation in large clinical and research cohorts is available in our preprint:

[**https://www.medrxiv.org/content/10.1101/2025.05.19.25327921**](https://www.medrxiv.org/content/10.1101/2025.05.19.25327921)

---

## **When to Use Talos**


Talos is designed to support **automated reanalysis of rare disease cohorts**, enabling identification of **candidate causative variants in known disease genes** based on the latest available evidence. It is best suited for scenarios where:

- You are performing **routine reanalysis** of undiagnosed individuals (e.g. monthly or quarterly)

- You want to detect **variants that have become reportable** due to updates in gene–disease or variant–disease knowledge

- You aim to **minimise the number of variants requiring manual review** optimising for specificity over sensitivity

- You are working with **exome or genome sequencing data** from previously analysed research or clinical cohorts

- You need a scalable, reproducible pipeline for **family-based or cohort-scale analysis**


Talos is **not currently designed** for:

- Identifying **novel candidate disease genes** or gene discovery

- Analysing **mitochondrial variants, short tandem repeats (STRs), mosaic variants**, or variants outside standard clinical reporting regions

> Support for some of these variant types may be added in future releases.

Talos complements existing variant curation workflows by focusing on high-specificity identification of variants that are likely to explain a participant’s condition, based on established gene–disease associations and up-to-date variant-level evidence.

---

## **🚀 Quick Start**

Talos is implemented using **Nextflow**, with all dependencies containerised via Docker. The example workflows can be run locally or on a cluster.

There are two primary workflows:

* `preparation.nf`: downloads and formats data in preparation for Talos runs
* `main.nf`: imports and executes the two sub-workflows which comprise the Talos runs

### **1. Install Requirements**

- [Nextflow](https://www.nextflow.io/docs/latest/install.html)

- Docker

To build the Docker image:

```
docker build -t talos:10.0.3 .
```

### **2. Download Annotation Resources**

Talos requires several large external resources (e.g. reference genome, gnomAD, AlphaMissense, Phenotype data). These are expected in a `large_files` directory. See [large_files/README.md](large_files/README.md) for detail on where to obtain them, and a [script](large_files/gather_file.sh) which will handle the initial download of all required resources.

### **3. Run Preparation Workflow**

In addition to the downloaded raw resources, Talos requires two other annotation sources to be kept up to date:

- ClinVar data, formatted into Hail Tables
- PanelApp data, an up-to-date dump as JSON

And a third data source (AlphaMissense) to be reformatted from the TSV into a Hail Table.

A separate sub-workflow, `preparation.nf` handles the download and formatting of this data:

```bash
nextflow \
    -c nextflow.config \
    run preparation.nf \
    [--processed_annotations <path>] \
    [--large_files <path>]
```

The parameter `processed_annotations` should point to a static directory where talos-generated files can be stored, and any future run of Talos will be able to access them. i.e. data prepared and written here is not linked to any individual underlying cohort or callset.

### **4. Run Annotation & Talos Combined Workflow**

> **NEW IN 10.0.0**
> Inputs for the Talos workflow are now provided in a single file, `--input_tsv`, instead of using several separate parameters.

The inputs for the Talos workflow are:
- **cohort**: a collective name to identify the input/results, used in output directory and file naming
- **path**: path to the Cohort's input data (VCF)
- **type**: type of the input data, see below
- **pedigree**: path to a Pedigree for the cohort, See details [here](docs/Pedigree.md)
- **config**: default available, path to the Talos config - see [example_config.toml](src/talos/example_config.toml) for an example, and the [Configuration README](docs/Configuration.md) for a full breakdown of all config parameters
- **history**: optional, path to previous results
- **ext_ids**: optional, path to ID mapping to present alternate IDs in the HTML report
- **seqr_map**: optional, path to ID mapping to generate hyperlinks to Seqr in the HTML report

The TSV file can contain any number of rows, each representing a distinct Cohort. A parallel Annotation & Talos run will be triggered for each input row, writing to a distinct output folder. An example TSV file has been provided to demonstrate.

The [annotation workflow](nextflow/annotation.nf) pre-processes and annotates variants. This workflow only needs to be run once per dataset, with the resulting MatrixTable(s) re-used with each iterative analysis.

#### Input Types 📂

The input TSV uses two columns to locate variant input; `path` and `type`. `path` is the location of the input file or directory. `type` is one of 3 values, **vcf, shards, ss_vcf_dir**

1. **vcf** a single multisample VCF. This will be split into shards and processed in parallel.
2. **shards** a directory of pre-sharded multisample VCF fragments, each shard containing all samples.
3. **ss_vcf_dir** single-sample VCFs, to be merged in the workflow, then sharded. These are detected using a glob, with the file extension controlled by `params.input_vcf_extension` (defaults to "vcf.bgz")

All results from the workflow will be written to a path pattern `{workflow.outputDir}/{cohort}_outputs`. This argument should point to a directory outside this repository, though for demonstration purposes the default is `./nextflow`.

The [main.nf](main.nf) workflow can be used to run both the main workflows, or where the data has been annotated previously, just the Talos workflow:

```bash
nextflow \
  -c nextflow.config \
  run main.nf \
  --input_tsv nextflow/inputs/test.tsv \
  -output-dir <path_to_output_dir>
```

>**For best results we advise repeating the Talos workflow on a regular cadence**
```bash
nextflow \
  -c nextflow.config \
  run talos_only.nf \
  --input_tsv nextflow/inputs/test.tsv \
  -output-dir <path_to_output_dir>
```

---

## **🔬 Input Validation**
The first step of the Talos workflow is a module called *StartupChecks*, which runs a number of input validations:

1. Checks a config file is present, and checks all required entries are present and have the correct type
2. Opens the Matrix Table and checks the schema and data types
3. Parses the Pedigree file and checks that it's well formatted and affected participants are present
4. Checks the ClinVar data, ensuring it is recent and has sufficient entries

This module will either run and complete, or run and fail, printing a collection of all encountered errors. If it fails, you will need to fix the errors before restarting the workflow.

---

## **⚙️ Configuration**

Talos as an application is configured through a single `TOML` file. This contains all thresholds and parameters for the steps of the Talos workflow. See [`example_config.toml`](src/talos/example_config.toml) as a baseline example, and [Configuration.md](docs/Configuration.md) for extended details on the role of each parameter, and its default value.

Talos as a Nextflow workflow is configured through configuration files, one each for the [annotation](nextflow/annotation.config) and [talos](nextflow/talos.config) stages of the workflow. These configurations define the cohort name, paths to input and annotation resources, and runtime parameters. [NextflowConfiguration.md](docs/NextflowConfiguration.md) contains a full description of the default values and role in the analysis.

## **📄 Outputs**

Talos produces structured outputs to support both manual review and downstream integration.

### **Primary Output:**

- *.json file listing all candidate variants for each proband

- Includes variant-level and gene-level evidence, inheritance checks, and phenotype match tags

### **Optional Outputs:**

- **HTML reports** summarising results for analysts or clinicians

- **Simplified TSV** for Seqr ingestion via MinimiseOutputForSeqr

### **Reanalysis Metadata:**

- first_seen: when the variant was first returned

- evidence_last_updated: when its evidence last changed

Only variants passing configured thresholds and logic modules are returned.

---

## **🔁 Reanalysis Mode**


Talos is designed to support **automated, iterative reanalysis** of undiagnosed cohorts. To do this it reads the results of previous analyses, and integrates them into the latest report. This is currently done by reading in prior analysis results, and incorporating the previous observations with each run. To use this behaviour, use the config setting `params.previous_results`. See [History](docs/Reanalysis.md) for more information.

### **How it works:**

1. Run full annotation + prioritisation once

2. In future cycles, keep ClinVar / PanelApp up to date **using the prep workflow**

3. Rerun prioritisation to return **newly supported variants**

By integrating the results of previous analyses with each new run, each variant in the output includes:

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

Each module can be configured through the `.toml` config file (see [Configuration.md](docs/Configuration.md))

---

## **📓 Citation**


If you use Talos in your research or clinical workflow, please cite:

> Welland MJ, Ahlquist KD, De Fazio P, et al. _Scalable automated reanalysis of genomic data in research and clinical rare disease cohorts._ medRxiv 2025.05.19.25327921; https://doi.org/10.1101/2025.05.19.25327921


BibTeX:

```
@article{welland2025talos,
  title     = {Scalable automated reanalysis of genomic data in research and clinical rare disease cohorts},
  author    = {Welland, Matthew J and Ahlquist, KD and De Fazio, Paul and Austin-Tse, Christina and Pais, Lynn and Wedd, Laura and Bryen, Samantha and Rius, Rocio and Franklin, Michael and Hall, Giles and et al.},
  journal   = {medRxiv},
  year      = {2025},
  doi       = {10.1101/2025.05.19.25327921},
  url       = {https://www.medrxiv.org/content/10.1101/2025.05.19.25327921},
}
```
