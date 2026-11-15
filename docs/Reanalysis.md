# Putting the "Re" in "ReAnalysis"

With Talos we have created a process which can both analyse data as a one-time event, and run re-analysis, building on previous results with each run. By referencing and building on previous results, we can prioritise new variants (novel classifications, or changed evidence) and reduce the amount of work needed to re-analyse a case by allowing users to filter out variants which were seen previously and have not changed.

Talos does this through storing a representation of each run's results, and feeding those forwards into future runs. By incorporating the variants and timestamps from a series of executions, we gradually build a record of all previously seen results, each with the original date of its observation. This allows us to see when a variant was first classified, and how its classification has changed over time. This is then reflected in the report, where the date shows the most recent observation of changed evidence.

For NextFlow, this is mediated through the `history` column in the input TSV.

* We advise that with each completed run, the `history` parameter is updated to the latest full_results JSON in the output folder.
* If the file does not exist or was not provided, all variant discovery dates will be set to the time of the current run.
* To make that update mechanical, every `main.nf` run writes `talos_input_YYYYMMDD.tsv` to the root of the output directory. This is a copy of the run's input TSV with each row's `history` cell repointed at the results JSON that run published for that cohort (`{cohort}_analysis_YYYYMMDD/{cohort}_full_results_YYYY-MM-DD.json`); every other column is carried over untouched, and a missing `history` column is appended. Review it, then use it as `--input_tsv` for the next cycle.

> NOTE! To clarify, Talos does not automatically perform re-analysis, and has no internal 'state'. For each run to iterate on the previous run's results, they must be fed into the new analysis using the input TSV. Providing previous result files will cause dates and variant observations to be carried across runs, without that each run is a 'first pass'.

## Process

During the run, Talos analyses data as standard. Once the result set has been generated, the previous file is loaded up and current results are compared to its contents:

* if a variant is newly detected, the dates are unchanged
* If a variant was seen before:
    * `first_tagged` is the latest of either:
      * the first date the variant was seen at any category
      * the first time the gene was Green in PanelApp
    * `date_of_phenotype_match` is None, if there is no phenotype match, otherwise it is set to the earliest date a phenotype match was observed
    * `evidence_last_updated` is the latest of:
      * the dates a category was assigned for the first time
      * the latest date an associated ClinVar record increased its number of stars
      * a variant was first a phenotype-match
      * a gene was newly green in PanelApp
