# PanelApp

PanelApp is a crowd-sourced gene curation knowledge base, aggregating information on gene-disease
associations from a wide range of specialists. Within PanelApp there are a number of panels, each
containing a set of genes associated with a disease, or a theme. Within a panel genes are rated using
a traffic light system:

* Red (minimal evidence)
* Orange (intermediate evidence)
* Green (strong evidence linking gene to theme)

For this project, we will focus analysis on a specific list of genes associated strongly with mendelian
disease, with an option to supplement this gene list with any number of additional PanelApp entries.
The Gene list we will use as the spine of the analysis is the [Mendeliome](https://panelapp.agha.umccr.org/panels/137/).

During an analysis run, the [Query Script](../reanalysis/query_panelapp.py) is used to pull down the latest data
for the Mendeliome panel (gene ENSG, symbol, and Mode of Inheritance). This is saved as a dictionary indexed on
GRCh38/Ensembl 90 ENSG.

## Additional Panels

If additional gene panel IDs are requested, the following steps are followed:

* the new panel's content is retrieved from PanelApp
* the Mendeliome (base) data is updated per-gene with new content
  * if a gene was on the mendeliome, add a `flag` with the additional panel ID/name
  * if a gene was on the mendeliome with no MOI, and the subsequent panel contains an MOI, update the MOI
  * if a gene was not on the mendeliome, add a new gene entry with the panel ID as a flag
* if multiple additional panels are added, this is repeated for each (extending `flags` for each panel a gene overlaps with)

## New Genes

Once all required panels are merged in, the final [optional] step is to check for 'new' genes. This can be done by
providing a gene list; a JSON file containing a single list of Strings. Each gene's entry will be annotated with
`new=True` if the PanelApp gene doesn't appear in the provided gene list, otherwise `new` remains at the default value
which is False.
