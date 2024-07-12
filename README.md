# Talos

[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff) ![test](https://github.com/populationgenomics/automated-interpretation-pipeline/actions/workflows/test.yaml/badge.svg) ![black](https://img.shields.io/badge/code%20style-black-000000.svg)

## Purpose

A variant prioritisation tool, aiming to assist clinical curators by sifting through large datasets and highlighting a
minimal number of plausibly disease causative genetic variants with high specificity.

The tool is designed to run repeatedly on large cohorts, so that gradual changes in functional annotations, gene-disease
associations, population frequencies, and variant calling pipeline improvements will cause new variants to surface over
time. As and when new variants of interest arise, curators can be alerted without having to manually repeat the cohort
analysis in full.

Talos aims to fast-track interpretation for the growing number of sequenced cohorts where curators lack time to complete
reanalysis work alongside the primary analysis for new clinical data. In light of numerous publications showing the
clinical utility of periodic reanalysis, this tool aims to lighten the workload on clinical analysts, while acting with
very high specificity to maximise the time efficiency of results interpretation.

## Strategy

Talos runs an analysis in two separate phases

1. Filter and categorise variants, identifying which deserve further processing based on consequence annotations.
2. Check each of those variants against the family structure of the participants in which it was found.

Variants only reach the final results both expected to be damaging, and the inheritance pattern shown fits the MOI
associated with the gene it is found in.

## Categories

The variant labelling stage of Talos implements a number of independent categories. Each category represents a decision
tree, using variant annotations to decide if a category label should be assigned. Each of these categories has been
designed to represent a block of curation logic - "if these criteria are all fulfilled, this could be relevant to
diagnosis"

These categories are each independent, providing a framework for adjustment, configuration, or extension to include more
variations in the future.

### Category 1

![CategoryBoolean1](design_docs/images/Category1.png)

If variant is rated as Pathogenic/Likely Pathogenic in ClinVar, minimum 1 'gold star' of associated evidence, we want to
flag that variant for review.

This category is exceptional in the sense that `Cat.1` variants are always processed under a partial-penetrance model -
even if the variant isn't a strict fit with the family or phenotype, we would want to be alerted (e.g. to look for a
second-hit by another method)

### Category 2

![CategoryBoolean2](design_docs/images/Category2.png)

A key reason for recovered diagnoses during reanalysis is the evolution of gene-disease understanding over time. This
Category aims to identify these variants by carrying state from the previous analysis (see below) and flagging where a
variant of at least moderate impact is newly associated with a disease. [`High Impact` consequence](
http://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html) here is based on the VEP definition of
`HIGH` consequence.

### Category 3

![CategoryBoolean3](design_docs/images/Category3.png)

This Category leverages the work of [LOFTEE](https://github.com/konradjk/loftee), a tool for identifying variants likely
to create loss of function with high confidence. When reviewing variants we require a `High Impact` consequence is
present, combined with either LOFTEE or a Clinvar P/LP rating (any number of stars).

### Category 4

![CategoryBoolean4](design_docs/images/Category4.png)

Here we accept a milder consequence (
any [`HIGH`consequence](http://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html) + `Missense`), but
only when present in the family with evidence of being _de novo_. This leverages the built-in [_de
novo_ method in Hail](
https://hail.is/docs/0.2/methods/genetics.html#hail.methods.de_novo), which is itself an implementation
of [Kaitlin Samocha's _de novo_ caller](https://github.com/ksamocha/de_novo_scripts). This implementation doesn't
naively accept the trio genotypes, but also applies some probability modelling to the genotype likelihoods, and searches
for alt-supporting reads at low levels in parents, before treating a _de novo_ variant call as validated.

### Category 5

![CategoryBoolean5](design_docs/images/Category5.png)

A simple category, here we use [SpliceAI](https://www.cell.com/cell/pdf/S0092-8674(18)31629-5.pdf) to identify
variants with a strong possibility of disrupting/adding splice junctions. This category will come under fire if the
underlying SpliceAI tool becomes monetised.

### Category 6

![CategoryBoolean6](design_docs/images/Category6.png)

A simple cateogry, here we use [AlphaMissense](https://www.science.org/doi/10.1126/science.adg7492) to identify
missense variants predicted to have a strong effect on the folded protein. A variant passes this test if the AM-assigned
'class' is `likely_pathogenic`. We have ambitions to make a second category here with a higher threshold applied to the
AM continuous score, instead of taking the ~0.56 threshold AlphaMissense natively uses to determine likely pathogenic.

### Category PM5

A little bit of secret sauce here - a piece of work twinned with Talos involved developing
[ClinvArbitration](https://github.com/populationgenomics/ClinvArbitration) - a re-summary of ClinVar data using altered
heuristics to aggregate multiple submissions for a variant. After creating new ClinVar results, we annotated the
pathogenic SNVs with VEP, and then re-index the results on `Protein & Codon`. In line with the ACMG evidence criteria
PM5 (this variant is a missense, and another missense at this same codon has a Pathogenic rating in ClinVar), we use
these
re-indexed results at runtime to apply the PM5 category.

This category is applied in the form `categorydetailsPM5=27037::1+27048::1`

- this is a `+` delimited list of entries (can be null)
- each entry is `ClinVar allele ID` :: `ClinVar Star rating`

This is processed upon variant ingestion, and back filtered to remove any associated ClinVar IDs which are this exact
variant.

## Reanalysis

The heart of Talos's utility is in enabling explicit re-analysis, which is done in two key ways:

### Gene Panel

* each time we query for gene panel/ROI data from [PanelApp](https://panelapp.agha.umccr.org/), we record the results
* if a gene features on a panel where it was previously absent, all variants in the gene will be eligible
  for `Category2` for the current run (if the participant in question has the panel applied)
* if a panel is applied in this analysis when it was previously not used, all variants in all genes on that panel are
  eligible for `Category2` for the current run (if the participant in question has the panel applied)

As we only search for 'Green' genes in PanelApp (those with sufficient evidence of disease association), a gene being
upgraded from Red or Amber-rated to Green will be picked up as a new gene.

Example:

* In the previous analysis, panel number 42 (phenotype: Boneitis) was applied, containing geneX and geneY
* In the current analysis, panel 42 newly features geneZ
* PatientA was phenotype-matched to panel 42, and has a geneZ variant marked as Category2, so this can reach the report
* PatientB was not phenotype-matched to panel 42, but has a geneZ variant marked as Category2. This category is stripped
  off, as panel 42 was not applied to this participant.

### Variant Results

* each time Talos runs, a minimised representation of the results is made; for each participant:
    - list the variants that were reported
    - the categories those variants were annotated with
    - and the date the category was first applied
    - the most recent date that the evidence changed (new category labels were applied)
* before the current report is written to file, the latest history file (if one exists) is checked:
    - if this variant & category was seen before, the `first_tagged` date is taken from the history file
    - if the variant was Pathogenic in ClinVar, the number of evidence stars is recorded. In future we may want to
      highlight improved evidence in the report as a review priority
    - if a category is assigned for the first time, the `evidence_last_updated` date is set to `today`

Together, these allow us to create reports which are easily filtered for events which occurred for the first time in
this latest run, removing all variants which were previously flagegd.

n.b. we do not currently hard-filter these results to remove previously-seen, we just enable that action by others

Example:

- If a variant appears as only `Cat.1`, and was previously a `Cat.1`, it will have `first_tagged`
  and `evidence_last_updated` set to the date of first appearance
- If a variant has been seen as a `Cat.1` before, and now is both `Cat.1` & `Cat.2`, `first_tagged` will be set to the
  date of first appearance, but `evidence_last_updated` will be set to `today`. The history file will be updated to show
  that `Cat.2` was applied `today`
- If a variant was never seen before and is now a `Cat.1`, `first_tagged` and `evidence_last_updated` will be set to
  today, and the history file will be updated to show that this variant was seen as a `Cat.1`, `today`
