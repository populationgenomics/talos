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

[8.3.5] - 2026-01-14

### Changed

* Talos now has a dependency on Mendelbrot, an external library which is acting as a common utilities collection for Talos and TalosAf
* Uses the PedigreeParser from Mendelbrot, deletes the duplicated code and tests in this codebase

### Added

* Adds HGVS interpretation to the amino_acid_change field coming out of BCFtools CSQ

[8.3.4] - 2025-11-27

### Fixed

* Substantially alters the de novo detection method, inferring GT (where missing but high quality -> HomRef), AD (where missing but high quality HomRef insert dummy values, where AD is a single element list add a 0 representing no alt reads), and DP (If DP is defined, use it, otherwise derive DP from AD)
* These changes should allow collaborators with gVCF combined data, but not combined with Hail Query, to run de novo detection

### Changed

* The docker images, both CPG internal and general purpose, both use BCFtools 1.22, and the NF CSQ annotation module updates to use appropriate syntax

[8.3.0] - 2025-10-18

### Fixed

* CPG Internally - due to altered CSQ behaviour in BCFtools 1.22, revert back to 1.21. Retain the explicit prefix unification in the CSQ call in case I accidentally update the version again.

[8.2.0] - 2025-09-22

### Added

* Various report reformats to improve usability and clarity.
* Results report now contains tabs for: “Results”, “Run Metadata”, and “Run Configuration”.

### Changed

* Variant table has been reformatted to support display on smaller screens.

[8.1.0] - 2025-09-19

### Added

* New default variant category:
  * `ClinVar Recent Gene` -  P/LP ClinVar with 0 stars that are in "new" PanelApp genes. New genes are defined by the recency threshold configured via GeneratePanelData.within_x_months (default: 24).
* New support/secondary category:
  * `ClinVar 0-star` - P/LP ClinVar with 0 stars. This is a support/secondary category by default, but can be promoted to strong in config.

### Changed

* GeneratePanelData.within_x_months is now 24 months by default (was 6 months). This means that "new" genes are those added to PanelApp within the last 24 months.

### Fixed

* Corrected boolean precedence in Hail expression for `categorybooleanclinvar0starnewgene`.

[8.0.4] - 2025-09-19

### Changed

* the Ensembl GTF file is updated to v115, this has been reflected in the downloader script and config.
* the list of transcript consequences where we are searching for de novo variants has been expanded to include in-frame deletion and insertion

[8.0.2] - 2025-09-03

### Removed

* The `result_history` mechanic for providing a directory of discrete files represnting state

### Changed

* State is now handled by passing a whole previous analysis result JSON into the workflow using `--previous_results`
* During MOI evaluation, a Homozygous variant will be considered for Comp-Het analysis, but only if the second hit is a SV deletion
* The ClinVar data download has been updated to the latest Zenodo release (2025-09)

[8.0.0] - 2025-08-19

### Removed

* All phenopacket generation and parsing code and dependencies are removed
* Removes the second report from HPOFlagging, instead exporting a single, phenotype-annotated report

### Changed

* No more phenopackets, instead we're using a single TSV (pedigree) file to contain the pedigree data, and optionally HPO terms. See the docs/Pedigree.md file for details on the new format and corresponding parser.
* The Report now presents the variant Categories as a short phrase (e.g. "ClinVarP/LP, HighImpact", instead of "1, 3")
* Lots of Model updates, but the liftover methods should accommodate them

### Added

* Added a Startup validation module to the Nextflow Talos-workflow. This checks the provided files/data should give successful results
* Added a pedigree parser module, and uses it internally
* Adds a file download process for all the required large input files
* Adds a breakdown of count/mean/median/SD for instances of each variant category in the report JSON

[7.5.0] - 2025-08-07

### Changed

* De Novo variant detection/filtering

During the de novo filtering process, GQ requirements have been altered. Instead of a flat "GQ >= 25" test against all
members of the callset, we now  have two separate configurable parameters:

* `de_novo.min_proband_gq` is filtered against all *affected* members of the pedigree.
* `de_novo.min_all_sample_gq` is optional, and is levied against every sample in the callset. Default is not to apply.

This altered behaviour means that we can now apply a stricter test against the proband, while still allowing
for a more relaxed test against other samples. This is useful in cases where the proband has a high GQ, but the parents
have a lower GQ, and we want to detect the de novo call.

A key situation to be aware of here - when combining single-sample VCFs into a multi-sample callset, as implemented in
the Nextflow workflow, WT calls with a GQ and AD of 0 are inserted when an individual had no evidence at the locus. This
means that if the `de_novo.min_all_sample_gq` test is applied, it will fail for all inserted WT calls, which will greatly
reduce candidate de novo discoverability. A logging message is printed to indicate when this test is run, and the config
example contains a comment indicating that this option will harm de novo discoverability if WT GQ=0 calls are inserted.

* Large files

Updated NextFlow config and README indicating that the `phenio.db` file is expected to be provided to the workflow
decompressed. A substantial portion of runtime and disk was spent decompressing this file, and keeping it decompressed
is a more efficient use of resources. The documentation/download process will have a substantial overhaul in an upcoming
change.

[7.5.0] - 2025-08-05

### Changed

1. The default NextFlow Glob for indidvidual VCFs has been changed from `*vcf.bgz` to `*.vcf.gz`, to match a more commonly used
   file extension The test data file names are updated to reflect the new default. This can be overridden with `--input_vcf_extension vcf.bgz`.
2. The NextFlow Glob also includes a `checkIfExists` check, so that if no VCFs are found, the workflow will not run and
   will instead throw an error. This is to prevent the workflow running with no input data, which was previously possible
   and would result in a confusing 0-second runtimes.
3. BCFtools merge now inserts WT genotypes for missing genotypes. This is to ensure the AC/AF/AN calculations are across
   the whole sample group under consideration. This is a change from the previous behaviour, which would not insert WT
   genotypes and would result in inflated AF values.
4. BCFtools no longer writes an uncompressed intermediate file after the merge. Instead it streams directly through the
   fill-tags plugin.
5. NF-workflow: Many stages have been changed from 'copy results to the output folder' to 'symlink to the output folder'.
   This is to reduce the amount of data copied around, and to ensure that the output folder is not bloated with various
   intermediate files which have no residual value after the workflow has concluded. This is a change from the previous
   behaviour, which would copy all results to the output folder. Depending on deployment, this may cause issues, and the
   copy behaviour may be preferable.
5. CPG-internal: no HT/MT outputs are tar'd to transmit between Stages; instead we use gcloud to copy the un'tar'd MT/HT.
   This action alone was adding ~35% to runtimes, and was not necessary for the workflow to run correctly.

### Removed

* The HpoFlaging step no longer writes two separate reports (HPO-matches annotated, and hard filtered to only HPO-matches).
  Instead it writes only a single output with HPO-matches annotated, and filtering can be applied in the HTML report.
* The CPG-Flow stage for squashing the MT into a Tarball has been removed.

### Added

* Nextflow has been added to the CPG-internal Dockerfile as standard

[7.0.8] - 2025-04-22

### Changed

* Default NextFlow workflow is updated to allow for a pre-existing phenopacket file to be used. Also adds 3 options:

1. if phenopackets are provided, use that
2. if a pedigree with hpo terms is provided, use that to generate a phenopacket file
3. if neither is provided, use the pedigree to generate an empty phenopacket file

* the zenodo link and filename for the gnomAD reference file has been updated (same file, with `.zip` removed)

* internal workflow allows for a specific clinvar path to be provided - now all projects have permission to query fewgenomes

### Deprecated

* remove the FindGeneSymbolMap file, not used any more
* remove the get_logger (named logger instance) method, replace with universal `loguru.logger`

[7.0.6] - 2025-04-10

### Fixed

CreateTalosHTML now correctly populated hyperlinks in the full-sized, and sample-subset, reports.

### Added

Nextflow annotation workflow now allows for users to provide a pre-merged VCF, skipping the merging step.

[7.0.5] - 2025-04-09

### Changed

The GeneratePanelData module now skips over any HPO terms which were not found in the ontology. This is an escape for
deprecated terms which are present in the metadata, previously these were causing a crash in the pipeline. This means
that some outdated HPO terms will not have a chance at a gene panel match, with logging to indicate which terms were
causing problems.

### Removed

setup.py has been deleted, in favour of a fully pyproject-based setup. This is a more modern approach to packaging, and
was in use as-of 7.0.0... I just forgot to remove the old setup.py file.

[7.0.4] - 2025-04-03

What happened to 7.0.3? Who knows!

### Changed

The way hyperlinks are generated out to external resources has been changed to be less Seqr-specific and more flexible.
Package is now limited to Python 3.10.X, due to a restriction in CPG-Flow.

Extensive renovation of how population filters are applied to different types of variant/MOI. This new approach
separates Dominant/non-Dominant filters, and ClinVar/Non-ClinVar filters. Each MOI model implents one of these 4 models,
each populating its filters from a specific block in config. As well as resulting in a more transparent and extensible
filtering process, this also removes a number of noisy variants we were retaining due to overly lax/absent filtering of
ClinVar Pathogenic variants.

### Added

A CPG-Flow workflow to run this pipeline internally at CPG. This is a wrapper process migrated from the CPG's internal
production-piplines workflow, with a few bits of updated syntax.




[7.0.2] - 2025-03-25

### Fixed

Correction to some of the NF modules - MatrixTables are no longer being re-compressed, so updating the expected file
extension to `.tar` instead of `.tar.gz` in the workflow.

[7.0.1] - 2025-03-24

### Fixed

Patched some SV VCF parsing - I was expecting the same schema to be present in CNV and GATK-SV VCFs, but some fields
need to be optional for both to pass through the parser.

### Added

A new block of functionality allows for genes and corresponding MOIs to be added from configuration, rather than only
pulling from PanelApp. This is useful for genes which are not in PanelApp, or for genes which are in PanelApp but not
expected to be found through the existing phenotype-matching process. It's also able to apply more lenient MOIs for
genes which were found in PanelApp, e.g. to trigger a more exploratory analysis.

[7.0.0] - 2025-03-20

### Origin

This is the starting point for Talos as an easily redistributable tool. Prior to this point Talos was reliant on
external annotations in a fixed format, making redeployment a challenge. This version is the first to own all annotation
processes, implementing an annotation pipeline which is minimal and quick to run.

The aim starting from this version is to ensure all Talos developments can be rolled out easily to external users with
no breaking changes, and no dependencies on external data formats or privately configured workflows.

### Added

Two nextflow workflows in the `nextflow` directory:

1. annotation.nf - a workflow to run the annotation pipeline, from raw VCF to annotated MT
2. talos.nf - a workflow to run the Talos pipeline, from annotated MT to report

These may not scale fantastically, but they are a good starting point for running Talos in a reproducible way, and they
can be used as a functional demo using test data distributed with the repository.

This annotation workflow is solely for SNV/Indels, but the SV analysis now runs from VCF instead of from the MT format
we hold internally. The annotations we expect are still in the exact format applied by GATK-SV's AnnotateVcf module, but
this is a good foundation for allowing SV VCFs annotated from other tools (e.g VEP to be used).

This annotation pipeline makes use of gnomAD V4.1 data, which is the latest version of gnomAD at the time of writing.
This is a far larger dataset than the previous gnomAD 2/3, and is expected to provide stronger filtering power.

[6.0.0] - 2024-10-21

### Added

* MakePhenopackets stage, which generates a Phenopacket file from Metamist

### Changed

* Renamed `GenerateSeqrFile` to `MinimiseOutputForSeqr`
* Renamed `vcf_to_mt` to `VcfToMt`
* The 'extended PED' file format is no longer used - instead we use a normal 6-col PED file, and a separate file with
  phenotypic data and external IDs for each participant in a GA4GH-compliant Cohort/Phenopacket format
  * An adapter remains in [talos/CPG](src/talos/ConvertPedToPhenopackets.py) to convert the old format to the new format

### Deprecated

* GeneratePED has been deleted

[5.0.0] - 2024-07-08

### Changed

* Throwing all the source code in a `/src` folder
* Principal scripts are renamed, and accessible as command-line entrypoints

1. GeneratePED
2. GeneratePanelData
3. QueryPanelapp
4. RunHailFiltering
5. RunHailFilteringSV
6. ValidateMOI
7. CreateTalosHTML

* extensive removing of the pre-commit exclusions, and consequent reformatting to pass new tests
* documentation updated

[4.0.0] - 2024-06-09

### Changed

- Name AIP -> Talos
- Pre-commit tooling changed again (linting with Ruff instead of Black)
- Removes use of QoB, and relies on 'local' Hail Query operations
- Drops Peddy in favour of Peds for a more minimal install
- Splits metamist querying into a separate script, with all other downstream operations using an extended PED format (6
  cols + ext IDs + HPO IDs)

MOI testing

- Uses the altered pedigree representation for MOI
- Each variant is interpreted relative to the immediate Trio (i.e. a quad is interpreted as 2 separate trios, not as a
  single variant/variants solving all presentations in a family)

### Added

- An adapter to take a VEP-annotated VCF and restructure in Hail to match the VEP-from-JSON annotated version we
  natively expect from the seqr loader pipeline

### Removed

- Everything to do with VEP... This does not contain a vep installation, and no longer determines how/if VEP is to be
  run. Except for ClinVar reanalysis - still in the process of dog-fooding that from ClinvArbitration instead of
  reimplementing here
- Removes even an example of analysis outcome registration from this repository, not even we use that functionality...

[3.2.2] - 2024-04-12

### Changed

Huge re-working of the pre-commit tooling, resulting in reduced overall line count

[3.0.0] - 2023-12-08

This was a MASSIVE refactor affecting almost all parts of the codebase. Instead of characterising
data in various places using type hinting, I've now introducted Pydantic models to formally represent
the data used in each part of the application (see reanalysis/models.py). This has the benefit of
making the code more readable, and also allows for the use of Pydantic's validation and parsing.

A key note is that all serialised/deserialised data is now done through a Pydantic validation layer,
forcing the data to conform to the model. This may require some objects (e.g. state/history files)
to be reformatted, as during this work I identified some inconsistencies in the data.

### Added

* Data models for every aspect of the analysis
* Pydantic validation layer for all serialised/deserialised data
* Specifically, all variants inherit from VariantCommon, and Structural and Small Variants each have
  a fundamentally different object representation, with a common interface but some internal methods
  specific to each variant type. This is extensible (e.g. further inheritance for STRs)

### Changed

* All minimial-case representations of Variants or Reportable events used in unit tests are now
  represented as full Pydantic models for the corresponding data type
* So many little things, this was a +3000/-2000 line change

[2.1.0] - 2023-11-20

This change involved the introduction of SV data, and the co-processing of SV and Small Variant
data in a unified analysis process. The SV data format expected is based on the GATK-SV pipeline
and annotation.

This also includes the addition of the `sv1` category, which represents LoF structural variants

### Added

* SV data is now processed alongside small variants
* `sv1` category specific to SV data

### Changed

* Some of the AbstractVariant loading is modified to allow for SV data (e.g. defaults for depths)

[2.0.0 & 2.0.1] - 2023-10-26

This bump is less substantial than the version number suggests, but in hindsight the previous
version should probably have been a major bump due to wholesale changes in how input is provided.

### Added

* CategoryBoolean6 (AlphaMissense) is now active, and should be implemented in a backwards
  compatible way (no failure if annotation is absent)
* VCF export includes `am_class` and `am_pathogenicity` (defaulting to empty String)
* If these fields are not missing, they should render in the report's variant drop-down drawer

### Changed

* The Hail Labelling stage also takes a dataset argument (used to generate a dataset-specific
  path for the temporary files/checkpoints)
* Removed a second config file from the test framework (consolidate all fields into 1 file)
* MyPy and Ruff are implemented as replacement linting tools. Stricter, and a billion times faster

### Deprecated

* Flake8 and PyLint are removed as linting tools

[1.2.0] - 2023-10-16

### Changed

* a lot. Like... a lot
* All the Cohort-specific content in the configuration files has been removed from the main config file. For CPG usage
  we have some sensitivities around the open publication of the cohorts and details of the analysis carried out. As a
  precaution (although this config does not contain sample-level data) we have removed the cohort-specific content from
  this repository and moved it to a private repository.
* Structurally the cohort-specific region of the config file has changed. This is to better facilitate running of Talos
  through our main analysis pipeline, requiring less manual intervention to run exomes and genomes separately. The
  utils methods `get_cohort_config` and `get_cohort_seq_type_conf` have been updated to reflect this.
    1. `config.cohorts.COHORT` contains details general to a whole cohort - blacklisted genes, additional panels, solved
       cases
    2. `config.cohorts.COHORT.GENOME` contains details specific to the genome analysis of a cohort - the Seqr project
       name, mapping between individual and seqr family IDs, and any labelled variants for this group
    3. `config.cohorts.COHORT.EXOME` is the counterpart to `.GENOME`
* Hail Labelling/Filtering now takes an output path, to be more easily included in a Hail Batch workflow
* The run preparation and HPO-panel script now use GQL queries, which drops a few methods
* The external_lookup file is no longer used - the same Internal-External ID mapping is noq stored alongside phenotype
  data in the HPO panel file.
* `--plink` argument to interpretation_runner is renamed to `--pedigree`
* Methods to reduce a cohort by a fixed percentage have been deleted - was useful as a potential test/train set, but has
  not been used recently

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
* See [this PR](https://github.com/populationgenomics/talos/pull/305)

[1.1.1] - 2023-05-08

### Added

* See [this PR](https://github.com/populationgenomics/talos/pull/293)
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
