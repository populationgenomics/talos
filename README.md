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
docker build -t talos:9.0.3 .
```

### **2. Download Annotation Resources**

Talos requires several large external resources (e.g. reference genome, gnomAD, AlphaMissense, Phenotype data). These are expected in a `large_files` directory. See [large_files/README.md](large_files/README.md) for detail on where to obtain them, and a [script](large_files/gather_file.sh) which will handle the initial download of all required resources.

### **3. Run Preparation Workflow**

In addition to the downloaded raw resources, Talos requires two other annotation sources to be kept up to date:

- ClinVar data, formatted into Hail Tables
- PanelApp data, an up-to-date dump as JSON

And a third data source (AlphaMissense) to be reformatted from the TSV into a Hail Table.

A separate sub-workflow, `talos_preparation.nf` handles the download and formatting of this data. Run this with:

```bash
nextflow -c nextflow.config run preparation.nf [--processed_annotations <path>] -with-report

 N E X T F L O W   ~  version 25.10.4

Launching `nextflow/preparation.nf` [irreverent_crick] DSL2 - revision: 6aacf78940

Attempting to download ClinVar raw data, requires internet connection.
If this step fails, try running gather_files.sh in the `large_files` directory.
executor >  local (2)
[b5/7767be] process > DownloadClinVarFiles (1)        [100%] 1 of 1 ✔
[f4/64b952] process > ResummariseRawSubmissions (1)   [100%] 1 of 1 ✔
[25/42f6eb] process > AnnotateClinvarWithBcftools (1) [100%] 1 of 1 ✔
[53/08b3b2] process > MakeClinvarbitrationPm5 (1)     [100%] 1 of 1 ✔
[f9/254c75] process > DownloadPanelApp (1)            [100%] 1 of 1 ✔
Completed at: 02-Feb-2026 14:50:20
Duration    : 5m 35s
CPU hours   : 0.2
Succeeded   : 3
```

The parameter `processed_annotations` should point to a static directory where talos-generated files can be stored, and any future run of Talos will be able to access them. i.e. data prepared and written here is not linked to any individual underlying cohort or callset.

### **4. Run Annotation & Talos Combined Workflow**

The [annotation workflow](nextflow/annotation.nf) pre-processes and annotates variants. This workflow only needs to be run once per dataset, with the resulting MatrixTable(s) re-used with each iterative analysis. The starting point for this workflow is one of the following:

1. a pre-sharded multisample VCF, with each shard containing all samples. Provide the directory with `--shards [path]`
2. a collection of single-sample VCFs, to be merged in the workflow. Provide the directory with `--ss_vcf_dir [path]`
3. a single multisample VCF. Provide this with `--vcf [path]`

The [talos workflow](nextflow/talos.nf) applies a series of filter- and labelling strategies to the annotated data, selecting variants likely to have clinical relevance.

All workflow outputs will be written to `--cohort_output_dir`. This argument should point to a directory outside this repository, though for demonstration purposes the default is `./nextflow/cohort_outputs`.

The [main.nf](main.nf) workflow can be used to run both the main workflows, or just the Talos workflow:

```
nextflow -c nextflow.config run main.nf \
  [--large_files <path>] \
  [--processed_annotations <path>] \
  [--vcf <path>] [--shards <path>] [--ss_vcf_dir <path>] \
  [--cohort_output_dir <path>] \
  --pedigree <path_to.ped>

```

```bash
nextflow -c nextflow.config run main.nf \
  --ss_vcf_dir nextflow/inputs/individual_vcfs \
  --cohort test \
  --pedigree nextflow/inputs/pedigree.ped \
  --container talos:9.0.3

 N E X T F L O W   ~  version 25.10.4

Launching `main.nf` [drunk_lumiere] DSL2 - revision: 65fae8c181

executor >  local (12)
[83/9bb229] process > ANNOTATION:MergeVcfsWithBcftools (1)       [100%] 1 of 1 ✔
[bd/82bc5a] process > ANNOTATION:SplitVcf (1)                    [100%] 1 of 1 ✔
[a3/392dff] process > ANNOTATION:NormaliseAndRegionFilterVcf (1) [100%] 1 of 1 ✔
[79/f47015] process > ANNOTATION:AnnotateWithEchtvar (1)         [100%] 1 of 1 ✔
[31/c4947e] process > ANNOTATION:AnnotateCsqWithBcftools (1)     [100%] 1 of 1 ✔
[3d/43c468] process > ANNOTATION:AnnotatedVcfIntoMatrixTable (1) [100%] 1 of 1 ✔
[96/cc5cff] process > TALOS:StartupChecks (1)                    [100%] 1 of 1 ✔
[0c/6e8396] process > TALOS:UnifiedPanelAppParser (1)            [100%] 1 of 1 ✔
[b8/710728] process > TALOS:RunHailFiltering (1)                 [100%] 1 of 1 ✔
[20/60f765] process > TALOS:ValidateMOI (1)                      [100%] 1 of 1 ✔
[ed/d594f4] process > TALOS:HPOFlagging (1)                      [100%] 1 of 1 ✔
[18/656d92] process > TALOS:CreateTalosHTML (1)                  [100%] 1 of 1 ✔
Completed at: 18-Feb-2026 16:13:08
Duration    : 2m 47s
CPU hours   : 0.1
Succeeded   : 12
```

For best results we advise repeating the Talos workflow on a regular cadence. In situations where the data has been annotated previously, and you only want to access the Talos sub-workflow:

```txt
nextflow -c nextflow.config run main.nf \
    --cohort test \
    --cohort_output_dir <path to annotated outputs folder> \
    --pedigree nextflow/inputs/pedigree.ped \
    -entry TALOS_ONLY

 N E X T F L O W   ~  version 25.10.4

Launching `main.nf` [ridiculous_shirley] DSL2 - revision: 65fae8c181

executor >  local (6)
[06/9a7fcd] process > TALOS_ONLY:TALOS:StartupChecks (1)         [100%] 1 of 1 ✔
[c8/8fe264] process > TALOS_ONLY:TALOS:UnifiedPanelAppParser (1) [100%] 1 of 1 ✔
[f0/f22e51] process > TALOS_ONLY:TALOS:RunHailFiltering (1)      [100%] 1 of 1 ✔
[c2/479946] process > TALOS_ONLY:TALOS:ValidateMOI (1)           [100%] 1 of 1 ✔
[de/64bd4e] process > TALOS_ONLY:TALOS:HPOFlagging (1)           [100%] 1 of 1 ✔
[6c/7b3220] process > TALOS_ONLY:TALOS:CreateTalosHTML (1)       [100%] 1 of 1 ✔
Completed at: 18-Feb-2026 16:58:39
Duration    : 1m 11s
CPU hours   : (a few seconds)
Succeeded   : 6
```

---

## **📥 Inputs**

Talos requires the following inputs:

| **Input type**         | **Description**                                                                                                                                                                                                                                                                                                                    |
|------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Variant data**       | Either: a set of individual sample VCFs, or a pre-merged multi-sample VCF. Using the example workflow, providing individual sample VCFs is only intended for relatively small numbers of samples per analysis. Using a  pre-merged, normalised multi-sample VCF will be  more efficient and scale to much larger sample numbers.   |
| **Pedigree file**      | A `tsv` file describing family structure, and optionally phenotypic terms per-participant. See details [here](docs/Pedigree.md).                                                                                                                                                                                                   |
| **Configuration file** | A `.toml` config file specifying all workflow settings. See [example_config.toml](src/talos/example_config.toml) for an example, and the [Configuration README](docs/Configuration.md) for a full breakdown of all config parameters                                                                                               |
| **Annotation files**   | Various, incuding [ClinvArbitration](https://github.com/populationgenomics/ClinvArbitration) data, gnomAD frequencies, Ensembl GTF, and MANE gene data. These are downloaded by running the [large_files setup script](large_files/gather_file.sh), or through the `talos_preparation` sub-workflow.                               |


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
