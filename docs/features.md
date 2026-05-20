# Features

Talos prioritises variants using a small set of rule-based **logic modules**, each aligned with specific ACMG/AMP evidence criteria. These modules are each present in code, implemented in Hail Query, and some can be modified by altering thresholds in configuration files. This page outlines the modules and the workflow features built on top of them.

---

## Variant logic modules

### Module types

- **Primary** modules are sufficient to trigger reporting on their own (e.g. `ClinVar P/LP`).
- **Supporting** modules are only used as second hits in recessive genes (e.g. `AlphaMissense`).

### Standard modules

| Module                            | Description                                                                                                                           |
|-----------------------------------|---------------------------------------------------------------------------------------------------------------------------------------|
| **ClinVar P/LP**                  | Pathogenic or Likely Pathogenic by ClinVArbitration.                                                                                  |
| **ClinVar Recent Gene**           | P/LP in a PanelApp "new" gene (became Green within the recency window configured by `GeneratePanelData.within_x_months`, default 24). |
| **High Impact**                   | Predicted high-impact protein consequences.                                                                                           |
| **De Novo**                       | Confirmed de novo in an affected individual.                                                                                          |
| **PM5**                           | Missense in a codon with a known pathogenic variant.                                                                                  |
| **LofSV**                         | Predicted loss-of-function structural variant.                                                                                        |
| **ClinVar 0-star** *(supporting)* | P/LP with 0 gold stars in ClinVar.                                                                                                    |
| **AlphaMissense** *(supporting)*  | AlphaMissense-predicted pathogenic missense variant.                                                                                  |

Each module is configurable through the `.toml` config file — see the [Configuration reference](Configuration.md).

!!! tip "Adding modules"
    The list above is the standard set; new categories can be added by following the [Adding New Categories](AddingNewCategories.md) guide.

---

## Reanalysis mode

Talos is designed to support **automated, iterative reanalysis** of undiagnosed cohorts. It reads the results of previous analyses and integrates them into the latest report. Enable this behaviour with the config setting `params.previous_results`. See [Reanalysis](Reanalysis.md) for full detail.

### How it works

1. Run full annotation + prioritisation once.
2. In future cycles, keep ClinVar / PanelApp up to date using the prep workflow.
3. Re-run prioritisation to surface **newly supported variants**, feeding in the previous run's results each time to identify the first date each sample~variant combination was observed.

By integrating prior results with each new run, every variant in the output includes:

- `first_seen` — original detection date.
- `evidence_last_updated` — last time the supporting evidence (ClinVar, PanelApp) changed.

> Reanalysis keeps manual review burden low by allowing users to filter to only those variants with newly actionable evidence in each cycle.

---

## Phenotype matching

Talos supports phenotype-driven filtering using **HPO terms**. See [Pedigree & Phenotypes](Pedigree.md) for how to provide phenotype data on the pedigree.

### Matching strategies

- **Patient-to-Gene** — semantic similarity between patient HPO terms and gene annotations.
- **Patient-to-Panel** — PanelApp panels are assigned to a patient when their terms match the panel's HPO tags.
- **Cohort-to-Panel** — manually assign panels to all individuals in a cohort via config.

When phenotype terms are provided, they are used to:

- Build a more accurate gene panel for each analysis by working with PanelApp to match disease-focused panels to HPO terms.
- Prioritise variants in the HTML report by highlighting variants in genes that are on disease-specific panels, or where the participant and gene HPO termsets share phenotypic similarities.

---

## Outputs

Talos produces structured outputs to support both manual review and downstream integration.

### Primary output

- A `*.json` file listing all candidate variants for each proband.
- Includes variant-level and gene-level evidence, inheritance checks, and phenotype match tags.
- One or more MatrixTables representing the annotated and reformatted input data. These are the starting point for subsequent runs of the workflow, skipping the annotation process.

### Optional outputs

- **HTML reports** summarising results for analysts or clinicians.
- **Simplified TSV** for Seqr ingestion via `MinimiseOutputForSeqr`.

### Reanalysis metadata

- `first_seen` — when the variant was first returned for the proband.
- `evidence_last_updated` — when supporting evidence last changed.

Only variants passing configured thresholds and logic modules are returned.

---

## Input validation

The first stage of every Talos run is a `StartupChecks` module that validates inputs before analysis begins. It either runs to completion or fails fast with a collected list of all encountered errors. See the [Getting Started](getting-started.md#what-happens-during-startup) page for the full list of checks.
