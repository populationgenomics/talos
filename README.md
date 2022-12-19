# Automated Interpretation Pipeline (AIP)

## Purpose

A variant prioritisation tool, aiming to assist clinical curators by sifting through
large datasets and highlighting a small number of plausibly disease causative genetic
variants with high specificity.

The tool is designed to run repeatedly on large cohorts, so that gradual changes in
functional annotations, gene-disease associations, population frequencies, and variant
calling pipeline improvements will cause new variants to surface over time. As and
when new variants of interest arise, curators can be alerted without having to
manually repeat the cohort analysis in full.

AIP assists clinical interpretation for the growing number of sequenced cohorts where
curators lack time to complete reanalysis work alongside the primary analysis for new
clinical data. In light of numerous publications showing the clinical utility of
periodic reanalysis, this tool aims to lighten the workload on clinical analysts,
while acting with very high specificity to maximise the time efficiency of results
interpretation.

## Strategy

AIP runs an analysis in two separate phases.
