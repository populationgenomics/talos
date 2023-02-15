# Clinvar Re-Processing

See [relevant development issue](https://github.com/populationgenomics/automated-interpretation-pipeline/issues/147)

## Context

ClinVar is valuable resource in identifying known Pathogenic & Benign variants within a
genomic dataset. By aggregating evidence from a range of submitters, we can utilise the
crowdsourced information to annotate current data with established clinical relevance.

ClinVar entries consist of:

* Individual Submissions, representing an assertion made by a submitter about the impact
of an individual allele.
* Allele summaries, which produce a top-line decision about each allele by aggregating all
relevant submissions.

During benchmarking of this application we have run into numerous instances of failing to
identify known pathogenic variants, due to conflicting ClinVar submissions results. On
closer inspection, the aggregation logic for ClinVar submissions seems too conservative
for our needs. An example of this is [here](https://ncbi.nlm.nih.gov/clinvar/variation/10/);
despite 24 Pathogenic submissions to only 2 Benign, the variant is given an overall status
of `Conflicting interpretations`. Whilst this is accurate, it obfuscates the bias towards
pathogenicity present in the individual submissions. When we annotate a dataset with ClinVar
consequences, all we have is this top-line decision, meaning that we are unable to flag such
variants for more manual scrutiny.

The role of AIP is not to make clinical decisions, but to identify variants of interest
for further review by analysts. In this setting we want to flag variants where manual
review of the submissions could signal a variant is worth consideration, even if it
doesn't appear pathogenic within the strict aggregation logic of ClinVar. To this end we
created a manual re-curation of ClinVar which:

* Allows for specific submitters to be removed from consideration (i.e. so that when we
run benchmarking analysis on cohorts, we can blind our ClinVar annotations to entries
originating from that cohort)
* Accepts an optional date threshold, to simulate a 'latest' summary at the given point in
time (discarding any submissions added or edited after the date, i.e. to simulate different
ClinVar time points using the same starting files)
* Defers to submissions after mainstream acceptance of ACMG criteria (estimated start 2016)
* Performs a more decisive summary, preferring a decision towards Pathogenic/Benign instead
of summarising any disagreements as `conflicting`

## Process

The re-summary is rapid, and can be repeated at regular intervals, taking the latest
available clinvar submissions each time it runs. The files used are the `submission_summary`
and `variant_summary` present on [the NCBI clinvar FTP site](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/).
We bin submissions into few discrete categories; see [the script (L37)](
../reanalysis/summarise_clinvar_entries.py) for the bins used.

1. Iterate over all individual submissions, removing any from blacklisted providers or
which were submitted after the date threshold. Collect all retained submissions per-allele.
2. For each allele, if any retained submissions were last edited after 2015 (representative ACMG
date), reduce submissions to only those. If no subs are from after 2015, retain all.
3. Find a summary 'rating' across all alleles, checking these scenarios until a match is found:

  * If an Expert Review/Clinical Guideline submission is present - choose that rating.
  * If both Pathogenic and Benign submissions are present, check for a confident majority (>= 60% in majority, <= 20% in minority). If there is a clear majority, choose that as the overall rating.
  * If both Pathogenic and Benign subs are present, but no clear majority, assign `Conflicting`.
  * If over half of submissions at the allele are `Uncertain`, rate as `Uncertain`.
  * If any Pathogenic submissions, take `Pathogenic`
  * If any Benign submissions, take `Benign`
  * No satisfied conditions - `Unknown/VUS`

4. A similar approach is followed for determining a `star` rating:

  * If any submissions are `Practice Guideline` -> `4 stars`
  * If any submissions are `Expert Review` -> `3 stars`
  * If any submissions `Criteria Provided` -> `1 stars`
  * `0 Stars`

At this stage we have each allele with a new summary and star rating. The allele ID is
matched up with the corresponding variant coordinates and ref/alt alleles from the
variant summary file, then the whole object is persisted as a Hail Table, indexed on
Locus and Alleles, ready to be used in annotation within Hail.