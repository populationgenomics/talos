# PanelApp

PanelApp is a crowd-sourced gene curation knowledge base, aggregating information on gene-disease
associations from a wide range of specialists. Within PanelApp there are a number of panels, each
containing a set of genes associated with a disease, or a theme. Within a panel genes are rated using
a traffic light system:

* Red (minimal evidence)
* Orange (intermediate evidence)
* Green (strong evidence linking gene to theme)

For this project, we will focus analysis on a specific list of genes associated strongly with mendelian
disease. The Gene list we will use is the [Mendeliome Gene Panel](https://panelapp.agha.umccr.org/panels/137/).

During an analysis run, the [Query Script](../reanalysis/query_panelapp.py) is used to pull down the latest data
for the mendeliome panel (gene ENSG, symbol, and Mode of Inheritance). This is saved as a dictionary indexed on
the ENSG

One key use case for this application is to chart the differences in gene list(s) over time. This script provides
two key ways to do this:

1. Provide a prior version as a command line argument. If this is done, the script will parse the Mendeliome's
content at the indicated version. Using this as comparison data, the 'latest' data is annotated with whether the
gene is newly green, or if the MOI has changed.
2. Provide a path to a gene list file. This should contain a single gene symbol per line, and the 'latest' data will
be annotated with 'new' if the PanelApp gene doesn't appear in the provided gene list
