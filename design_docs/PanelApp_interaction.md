# PanelApp

PanelApp is a crowd-sourced gene curation knowledge base, aggregating information on gene-disease associations from a
wide range of specialists. Within PanelApp there are a number of panels, each containing a set of genes associated with
a disease, or a theme. Within a panel genes are rated using a traffic light system:

* Red (minimal evidence)
* Orange (intermediate evidence)
* Green (strong evidence linking gene to theme)

For this project, we will focus analysis on a specific list of genes associated strongly with mendelian disease, with an
option to supplement this gene list with any number of additional PanelApp entries. The Gene list we will use as the
core of the analysis deafults to the [Mendeliome](https://panelapp.agha.umccr.org/panels/137/), though this can be
altered through the configuration file.

During an analysis run, the [Query Script](../reanalysis/query_panelapp.py) is used to pull down the latest core panel
data (gene ENSG, symbol, and Mode of Inheritance), saved as a dictionary indexed on GRCh38/Ensembl 90 ENSG.

## Per-participant Panels

This codebase contains a [HPO-Panel Helper](../helpers/hpo_panel_matching.py), integrating participant phenotypic data,
PanelApp per-panel relevant phenotypes, and the HPO heirarchy to find relevant panels to apply to each individual
participant. This will be difficult to replicate outside of CPG infrastructure, so this describes the output generated,
which can be produced by other means:

```python
# participantID here must match the ID in the joint call
_panel_dict = {
    'participantID': {
        'panels': [1, 23, 456]
    }
}
```

If this file is supplied at runtime, any panel matched to a participant will be integrated into the overall panel query.
Algorithmically, the following process occurs:

1 Query for the core panel (Mendeliome, unless altered in config)

  * We collect the relevant MOI and symbol for each Green gene

2 Gather all phenotype-matched panels from the file
3 Query for each extra panel in turn, adding new genes from the additional panels

  * At each stage, if MOI was missing in a prior panel and present in the current panel, update the MOI
  * Record all the panels each gene appears on

## Notion of `New`

As mentioned in the [README](../README.md#gene-panelroi), we can supply prior data to AIP to ensure reanalysis is as
specific as possible. This data consists of all previously seen genes, and all panels on which each gene has been seen.
A gene is considered to be 'new', when either:

* it is seen on the current run, and has never been seen before
* it is on a panel where it has not previously occurred

If either of these conditions are satisfied, Panel IDs where the gene is `new` are added to a list. The cumulative panel
data is updated to record all panels in which the current genes have previously been seen.

### Application of `New`

When a variant is provisionally classified as `Cat.2` (new and impactful), we check that the gene was new in a panel
which applies to that specific participant. If a participant has a variant in a gene `new` in the analysis, this is only
relevant if one or more panels where the gene is found to be `new` applies to this participant.
