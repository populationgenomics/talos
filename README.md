# Automated Interpretation Pipeline (AIP)

[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff) ![test](https://github.com/populationgenomics/automated-interpretation-pipeline/actions/workflows/test.yaml/badge.svg) ![black](https://img.shields.io/badge/code%20style-black-000000.svg)

## Purpose

A variant prioritisation tool, aiming to assist clinical curators by sifting through large datasets and highlighting a
minimal number of plausibly disease causative genetic variants with high specificity.

The tool is designed to run repeatedly on large cohorts, so that gradual changes in functional annotations, gene-disease
associations, population frequencies, and variant calling pipeline improvements will cause new variants to surface over
time. As and when new variants of interest arise, curators can be alerted without having to manually repeat the cohort
analysis in full.

AIP aims to fast-track interpretation for the growing number of sequenced cohorts where curators lack time to complete
reanalysis work alongside the primary analysis for new clinical data. In light of numerous publications showing the
clinical utility of periodic reanalysis, this tool aims to lighten the workload on clinical analysts, while acting with
very high specificity to maximise the time efficiency of results interpretation.

## Strategy

AIP runs an analysis in two separate phases

1. Filter and categorise variants, identifying which deserve further processing based on consequence annotations.
2. Check each of those variants against the family structure of the participants in which it was found.

Variants only reach the final results both expected to be damaging, and the inheritance pattern shown fits the MOI
associated with the gene it is found in.

## Categories

AIP uses the concept of a 'Category' to label variants which are anticipated to cause disease. These categories are
assigned on the basis of consequence annotations, and each has been defined in collaboration with Clinical Curators to
codify a mental model for when a variant is likely to be relevant to a diagnosis.

These categories are each independent, providing a framework for adjustment, configuration, or extension to include more
variations in the future.

### Category 1

![CategoryBoolean1](design_docs/images/Category1.png)

If variant is rated as Pathogenic/Likely Pathogenic in ClinVar, minimum 1 'gold star' of associated evidence, we want to
flag that variant for review.

This category is exceptional in the sense that `Cat.1` variants are exempt from MOI tests - we want to report these
variants even if the inheritance pattern doesn't fully fit within the family.

### Category 2

![CategoryBoolean2](design_docs/images/Category2.png)

A key reason for recovered diagnoses during reanalysis is the evolution of gene-disease understanding over time. This
Category aims to identify these variants, by carrying state from the previous analysis (see below) and flagging where a
variant of at least moderate impact is newly associated with a disease. [`High Impact` consequence](
http://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html) here is based on the VEP definition of
`HIGH` consequence

### Category 3

![CategoryBoolean3](design_docs/images/Category3.png)

This Category leverages the work of [LOFTEE](https://github.com/konradjk/loftee), a tool for identifying variants likely
to create loss of function with high confidence. When reviewing variants we require a `High Impact` consequence is
present, combined with either LOFTEE or a Clinvar P/LP rating.

### Category 4

![CategoryBoolean4](design_docs/images/Category4.png)

This Category is the only definition that includes an MOI requirement - here we accept a milder consequence (any [`HIGH`
consequence](http://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html) + `Missense`), but only when
present in the family with evidence of being _de novo_. This leverages the built-in [_de novo_ method in Hail](
https://hail.is/docs/0.2/methods/genetics.html#hail.methods.de_novo), which is itself an implementation of [Kaitlin
Samocha's _de novo_ caller](https://github.com/ksamocha/de_novo_scripts). This implementation doesn't naively accept the
trio genotypes, but also applies some probability modelling to the genotype likelihoods, and searches for alt-supporting
reads at low levels in parents, before treating a _de novo_ variant call as validated.

### Category 5

![CategoryBoolean5](design_docs/images/Category5.png)

Another simple cateogry, here we use [SpliceAI](https://www.cell.com/cell/pdf/S0092-8674(18)31629-5.pdf) to identify
variants with a strong possibility of disrupting/adding splice junctions. This category will come under fire if the
underlying SpliceAI tool becomes monetised.

### Category Support

![CategoryBooleanSupport](design_docs/images/CategorySupport.png)

A purely _in silico_ category, here we accept 'consensus' across a number of informative prediction tools, with the
exact thresholds influenced by published work by [Pejaver, _et al_](
https://www.biorxiv.org/content/10.1101/2022.03.17.484479v1). Unlike the other categories, this is not enough to mark
variants as worth investigation on its own - due to high False Positive rate, this category is only relevant as a
second-hit level of confidence (i.e. a biallelic MOI can be satisfied by a Cat.Support variant combined with a variant
of any other category)

## Reanalysis

The heart of AIP's utility is in enabling explicit re-analysis, which is done in two key ways. Each of these is enabled
through accumulating data over a number of runs, then leveraging cumulative data to remove redundant results. This helps
ensure a curator's time is used as effectively as possible, by not presenting the same variants each time AIP runs.

### Gene Panel/ROI

[PanelApp](https://panelapp.agha.umccr.org/) is the primary source of Gene information used in AIP. For each analysis we
query for the current content of the [Mendeliome Panel](https://panelapp.agha.umccr.org/panels/137/), containing all
genes associated with Mendelian disease, other than those associated with unintended diagnoses (see [the
complementary Incidentalome](https://panelapp.agha.umccr.org/panels/126/)). From this Panel, we obtain all 'Green' genes
and their associated Modes Of Inheritance (MOI). Optionally we can add additional panels which will be integrated
alongside the Mendeliome content.

Prior data can be provided in JSON format, containing all genes which have previously formed part of the ROI, and a list
of all panels on which they have previously been seen. A gene will be treated as `New` if either it has never been seen
before, or features in a phenotype-matched panel on which it has not been seen before. i.e.

- if a Gene is promoted to `Green` on the Mendeliome, it will be recorded as `New`, and the prior data will be extended
  to show that the gene was seen on the Mendeliome.
- if a Gene has previously appeared on the Mendeliome, and during current run it now appears on a new panel, the gene
  will be recorded as `New`, and the new panel will be added to the prior data list.

### Variant Results

Once the core process of AIP has completed, and the reportable results are assembled, the final (optional) stage is to
back-filter against previously seen results. This prior data is assembled on a per-cohort basis, and contains results
previously seen, indexed by Sample ID.

This data consists of the variant IDs, Categories they've previously been assigned, and any supporting variants they
have formed a compound-het with. When reviewing the variants of a current run, we check for previously seen variants on
a per-sample basis, e.g.:

- If a variant has been seen as a `Cat.1` before, and appears again as a `Cat.1`, it will be removed
- If a variant has been seen as a `Cat.1` before, and now is both `Cat.1` & `Cat.2`, the `Cat.1` assignment will be
  removed, and will be reported only as a `Cat.2`. The prior data will be extended to show that it has been
  a `Cat.1 & 2`
- If a variant was never seen before, it will appear on the report with no removed Categories
- If a variant was seen in a compound-het now, and was previously partnered with a different variant, all `Categories`
  will be retained, and the new partner ID will be added to the list of `support_vars` in the prior data
