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

[1.0.1] - as yet unreleased

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
