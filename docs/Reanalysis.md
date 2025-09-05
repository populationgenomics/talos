# Putting the "Re" in "ReAnalysis"

With Talos we have created a process which can both analyse data as a one-time event, and run re-analysis, building on previous results with each run. By referencing and building on previous results, we can prioritise new variants (novel classifications, or changed evidence) and reduce the amount of work needed to re-analyse a case by allowing users to filter out variants which were seen previously and have not changed.

Talos does this through storing a representation of each run's results, and feeding those forwards into future runs. By incorporating the variants and timestamps from a series of executions, we gradually build a record of all previously seen results, each with the original date of its observation. This allows us to see when a variant was first classified, and how its classification has changed over time. This is then reflected in the report, where the date shows the most recent observation of changed evidence.

For NextFlow, this is mediated through the `params.previous_results` config entry, or `--previous_results <result_JSON>`, which points to a result from a previous run. If the file does not exist or was not provided, all variant discovery dates will be set to the time of the current run. If a previous result is passed in, the current run's results will be updated to reflect the first time events were discovered, and which events are discovered for the first time in the current run.

## Process

During the run, Talos analyses data as standard. Once the result set has been generated, the previous file is loaded up and current results are compared to its contents:

* if a variant is newly detected, the dates are maintained as the current date
* If a variant was seen before:
  * `first_tagged` is set to the date of its first observation, of any category
  * `evidence_last_updated` is set to the most recent date a category was assigned for the first time
  * if a variant was seen before, and is now seen with a novel comp-het partner, `evidence_last_updated` is today
  * `date_of_phenotype_match` is None, if there is no phenotype match, otherwise it is set to the earliest date a phenotype match was observed
