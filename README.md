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

- You aim to **minimise the number of variants requiring manual review** without compromising sensitivity

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


### **1. Install Requirements**

- [Nextflow](https://www.nextflow.io/docs/latest/install.html)

- Docker

To build the Docker image:

```
docker build -t talos:7.4.1 .
```

### **2. Download Annotation Resources**

Talos requires several large external resources (e.g. reference genome, gnomAD, AlphaMissense, Phenotype data). These are expected in a `large_files` directory. See [large_files/README.md](large_files/README.md) for detail on where to obtain them.

### **3. Run Annotation Workflow**

This step pre-processes and annotates variants. This workflow only needs to be run once per dataset:

```
nextflow -c nextflow/annotation.config \
  run nextflow/annotation.nf \
  [--large_files <path>] \
  [--processed_annotations <path>] \
  [--merged_vcf <path>]
```

### **4. Run Variant Prioritisation**

After annotation is complete, run the Talos prioritisation workflow. This workflow is the step required to be re-run regularly for an updated analysis:

```
nextflow -c nextflow/talos.config \
  run nextflow/talos.nf \
  --matrix_table nextflow/cohort_outputs/cohort.mt
```

---

## **📥 Inputs**

Talos requires the following inputs:

| **Input type**                  | **Description**                                                                                                                                                                                                                                                                                                                  |
| ------------------------------- |----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Variant data**                | Either: a set of individual sample VCFs, or a pre-merged multi-sample VCF. Using the example workflow, providing individual sample VCFs is only intended for relatively small numbers of samples per analysis. Using a  pre-merged, normalised multi-sample VCF will be  more efficient and scale to much larger sample numbers. |
| **Pedigree file**               | A `.ped` file describing family structure ([specification](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format))                                                                                                                                                                                  |
| **Phenotype data** _(optional)_ | Phenotypes encoded as HPO terms in a [GA4GH Cohort](https://phenopacket-schema.readthedocs.io/en/latest/cohort.html). An example is provided [here](nextflow/inputs/test_cohort_phenopackets.json), and an extended [README here](docs/Phenopackets.md)                                                                          |
| **ClinVar data**                | Pre-processed using [ClinvArbitration](https://github.com/populationgenomics/ClinvArbitration). The most recent version is available on [Zenodo](https://zenodo.org/records/15896156)                                                                                                                                            |
| **Configuration file**          | A `.toml` config file specifying cohort name, paths, and filtering settings. See [example_config.toml](src/talos/example_config.toml) for an example, and the [Configuration README](docs/Configuration.md) for a full breakdown of all config parameters                                                                        |
|                                 |                                                                                                                                                                                                                                                                                                                                  |

> A stub ClinVar dataset is included for testing: [clinvarbitration.tar.gz](nextflow/inputs/clinvarbitration.tar.gz). This file contains only 4 variants and is **not suitable for real analyses**. For real-data runs please visit the above Zenodo link and download the latest complete dataset.

---

## **⚙️ Configuration**

Talos as an application is configured through a single `TOML` file. This contains all thresholds and parameters for the steps of the Talos workflow. See [`example_config.toml`](src/talos/example_config.toml) as a baseline example, and [Configuration.md](docs/Configuration.md) for extended details on the role of each parameter, and its default value.

Talos as a Nextflow workflow is configured through configuration files, one each for the [annotation](nextflow/annotation.config) and [talos](nextflow/talos.config) stages of the workflow. These configurations define the cohort name, paths to input and annotation resources, and runtime parameters. [NextflowConfiguration.md](docs/NextflowConfiguration.md) contains a full description of the default values and role in the analysis.

## **📄 Outputs**

Talos produces structured outputs to support both manual review and downstream integration.

### **Primary Output (per family):**

- *.json file listing candidate variants

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


Talos is designed to support **automated, iterative reanalysis** of undiagnosed cohorts.

### **How it works:**

1. Run full annotation + prioritisation once

2. In future cycles, update ClinVar / PanelApp only

3. Rerun prioritisation to return **newly supported variants**

By integrating the results of previous analyses with each new run, each variant in the output includes:

- first_seen: original detection date

- evidence_last_updated: last evidence update (ClinVar, PanelApp)


> Talos maintains low review burden by returning only variants with newly actionable evidence.

---

## **🧬 Phenotype Matching**



Talos supports phenotype-driven filtering using **HPO terms**.



### **Matching Strategies:**

- **Patient-to-Gene**: semantic similarity between HPO terms and gene annotations

- **Patient-to-Panel**: PanelApp panels assigned if patient terms match panel HPO tags

- **Cohort-to-Panel**: manually assign panels to all individuals in config


Phenotype data is optional, and can be provided as a [GA4GH Cohort](https://phenopacket-schema.readthedocs.io/en/latest/cohort.html) ([README here](docs/Phenopackets.md)). When enabled, phenotype matching is used to:

* Build a more accurate gene panel for each analysis, working with PanelApp to match disease-focused panels to HPO terms
* Prioritises variants in the report by visually indicating variants in genes which are on disease-specific panels, or where the participant and Gene HPO termsets share phenotypic similarities

---

## **🧠 Variant Logic Modules**


Talos prioritises variants using rule-based **logic modules**, each aligned with specific ACMG/AMP evidence criteria.


### **Module Types**

- **Primary**: sufficient to trigger reporting on their own (e.g. ClinVar_PLP)

- **Supporting**: used only as second hits in recessive genes (e.g. In_silico)


### **Standard Modules**

| **Module**  | **Description**                                     |
|-------------|-----------------------------------------------------|
| ClinVar_PLP | Pathogenic or Likely Pathogenic by ClinvArbitration |
| High_Impact | Predicted high-impact protein consequences          |
| De_novo     | Confirmed de novo in affected individual            |
| PM5         | Missense in codon with known pathogenic variant     |
| SV_LOF      | Predicted loss-of-function structural variant       |
| In_silico   | AlphaMissense-predicted pathogenic missense variant |

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
