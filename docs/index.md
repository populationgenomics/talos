# Talos

**Talos** is a scalable, open-source variant prioritisation tool designed to support automated reanalysis of genomic data in rare disease. It identifies **candidate causative variants in known disease genes** by integrating static annotations (e.g. population frequency, predicted consequence) with dynamic knowledge sources such as ClinVar and PanelApp Australia. Talos applies a set of configurable, rule-based logic modules aligned with ACMG/AMP criteria and prioritises variants consistent with expected mode of inheritance and, optionally, patient phenotype.

While Talos can be used for one-off reanalysis of individual families or cohorts, its core design is optimised for **routine, cohort-scale reanalysis**. By comparing current annotations with prior results, Talos highlights variants that have become reportable due to newly available evidence — such as new gene–disease or variant–disease relationships — since the last analysis cycle.

A full description of the method and its validation in large clinical and research cohorts is available in our preprint: [medRxiv 2025.05.19.25327921](https://www.medrxiv.org/content/10.1101/2025.05.19.25327921).

> Note - whether you are a new Talos user, or have an implementation already, we encourage you to run the `preparation` workflow with each update!

---

## Latest Changes

{%
   include-markdown "../CHANGELOG.md"
   start="<!--latest-start-->"
   end="<!--latest-end-->"
   comments=false
%}

[Full changelog](changelog.md)

---

## Where to next

<div class="grid cards" markdown>

- :material-rocket-launch: **[Getting Started](getting-started.md)**
  Install the requirements, download annotation resources, and run your first cohort.

- :material-puzzle: **[Features](features.md)**
  Logic modules, reanalysis mode, phenotype matching, and what Talos is (and isn't) for.

- :material-cog: **[Configuration](Configuration.md)**
  Full reference for the Talos TOML config and Nextflow parameters.

- :material-history: **[Changelog](changelog.md)**
  Release history and version-by-version changes.

</div>

---

## When to use Talos

Talos is best suited for scenarios where:

- You are performing **routine reanalysis** of undiagnosed individuals (e.g. monthly or quarterly).
- You want to detect **variants that have become reportable** due to updates in gene–disease or variant–disease knowledge.
- You aim to **minimise the number of variants requiring manual review**, optimising for specificity over sensitivity.
- You are working with **exome or genome sequencing data** from previously analysed research or clinical cohorts.
- You need a scalable, reproducible pipeline for **family-based or cohort-scale analysis**.

Talos is **not currently designed** for:

- Identifying **novel candidate disease genes** or gene discovery.
- Analysing **short tandem repeats (STRs), mosaic variants**, or variants outside standard clinical reporting regions.

> Support for some of these variant types may be added in future releases.

---

## Citation

If you use Talos in your research or clinical workflow, please cite:

> Welland MJ, Ahlquist KD, De Fazio P, et al. *Scalable automated reanalysis of genomic data in research and clinical rare disease cohorts.* medRxiv 2025.05.19.25327921; <https://doi.org/10.1101/2025.05.19.25327921>

```bibtex
@article{welland2025talos,
  title     = {Scalable automated reanalysis of genomic data in research and clinical rare disease cohorts},
  author    = {Welland, Matthew J and Ahlquist, KD and De Fazio, Paul and Austin-Tse, Christina and Pais, Lynn and Wedd, Laura and Bryen, Samantha and Rius, Rocio and Franklin, Michael and Hall, Giles and et al.},
  journal   = {medRxiv},
  year      = {2025},
  doi       = {10.1101/2025.05.19.25327921},
  url       = {https://www.medrxiv.org/content/10.1101/2025.05.19.25327921},
}
```
