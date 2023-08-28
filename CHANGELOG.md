# Changelog

All notable changes to this project will be documented in this file.

Suggested headings per release (as appropriate) are:

* `Added` for new features.
* `Changed` for changes in existing functionality.
* `Deprecated` for soon-to-be removed features.
* `Removed` for now removed features.
* `Fixed` for any bug fixes.
* `Security` in case of vulnerabilities.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

[1.1.5] - 2023-08-28

### Changed

* reanalysis/interpretation_runner.py creates two reports; the first report is all standard content from the run, the
  second report is variants new within this latest run
* Relevant interfaces with the output file registration and report hunter scripts are altered
* The report hunter script is altered to create two separate index pages; one for the standard reports, one to only look
  for the latest-variant report

[1.1.4] - 2023-06-28

### Changed

* Updated to use the latest models in Metamist/sample metadata DB
* Queries updated to use GraphQL

[1.1.3] - 2023-06-22

### Added

* Adds a helper script which exports a file formatted specifically for ingestion into Seqr

[1.1.2] - 2023-05-31

### Changed

* Adds a number of changes to how ClinVar pathogenic variants are processed
* See [this PR](https://github.com/populationgenomics/automated-interpretation-pipeline/pull/305)

[1.1.1] - 2023-05-08

### Added

* See [this PR](https://github.com/populationgenomics/automated-interpretation-pipeline/pull/293)
* Adds the PM5 category, and the VEP content to create the relevant annotation data

[1.0.1] - 2023-04-01

### Added

* Added this changelog
* Allow for addition of non-Metamist metadata registrars, selectable in config
* Added Variant Tagging in report, adopted from MW@MSFT's fork - #271
* Add family ID into top level variant table

### Changed

* Separated out metadata registration into a separate file & job
* File registration is set to `always_run`, if indicated in config
* De Novo calling (hail method) can have `max_parent_ab` ratio overridden in config
* Establishes a method-contract for altering metadata registrars and HTML scripts
