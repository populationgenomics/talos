# HPO and Panel Matching

As part of the move towards personalisation of reported variants, we are incorporating phenotype information for the
analysed cohorts. The steps of this process are:

1. Ingesting the participant data, including the HPO terms associated with each individual
2. Querying PanelApp for all panels which have associated HPO terms
3. Digest HPO tree file [current source is hp.obo](https://hpo.jax.org/app/download/ontology) into a query-able format
4. Match between cohort HPO terms and HPO terms in PanelApp to find a group of panels relevant to each individual
5. Exporting this data in a format which can be used in AIP


## PanelApp Querying

We are querying the overview endpoint of the PanelApp API, which provides brief details on each panel, to gather
metadata. This endpoint doesn't contain the full gene list content of each panel, but it does express the 'related
panels' information which is being populated with HPO data by curators.

For each panel, we are using a regex match to find any valid HPO terms, and storing each panel ID linked to the relevant
HPO terms (so that one HPO term can be linked to multiple panel IDs without redundancy/duplication). This mapping is
then returned as a dictionary.

Note - there is an HPO field in the PanelApp database, but that data is not currently being exported (see [this issue](
https://gitlab.com/genomicsengland/panelapp/panelapp/-/issues/24))

## HPO Ontology parsing

The HPO ontology tree is downloaded in `.obo` format, which can be digested as a graph. In this format, each node can be
queried directly, as well as being linked to all parents/children. [obonet](https://pypi.org/project/obonet/) is used to
ingest the ontology into a [networkx](https://networkx.org/) graph, where all data can be queried.

## Participant Searching

Currently we are geared up to pull data from metamist, the CPG metadata API. The client for this is currently broken, so
we are using a dump instead.

From this dump, we are finding ParticipantIDs, FamilyIDs, and any corresponding HPO terms.

## HPO to panel searching

Get all unique HPO terms across the cohort, then for each unique HPO term:

1. If the HPO term has an exact match from the PanelApp data; if so, add to a set
2. If the term has a parent (less specific) recursively call the search using the less-specific term
   - If a term is deprecated or obsolete, check for any replacement terms listed in the ontology and re-run the search
   - if a term has 'alternative' terms, re-run using those terms
3. Repeat this search, becoming increasingly vague, until reaching a configured limit
4. Once the limit is reached, return the collection of all matched panel IDs

At the end of this process, have a lookup of all HPO terms found in the cohort, matched to all panels having HPO terms
within `n` layers of the HPO ontology. Alternatively we could continue all the way up to the root term, assuming that
terms associated with individual panels are specific.

## Matching panels back to participants

The final stage is to associate these identified panels back to the participants. This follows the following algorithm:

1. For each participant in turn, pull out the list of HPO terms
2. For each HPO term, find any panels which were fuzzy-matched to that term
3. Store the associated panels alongside the HPO term data for the participant


## Result

The end product is a dictionary indexed primarily on Participant ID, with each participant linked to their associated
HPO terms, panels to use, and Family ID.
