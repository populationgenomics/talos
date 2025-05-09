# Talos

[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff) ![test](https://github.com/populationgenomics/automated-interpretation-pipeline/actions/workflows/test.yaml/badge.svg)

## Overview

Talos is a Python variant prioritisation tool which finds diagnostically relevant variants within callsets.

It incorporates consequence annotation, family structures, participant phenotypes, and crowdsourced clinical knowledge
to identify variants likely to cause participant phenotypes. It does so with high sensitivity, whilst retaining maximal
specificity to reduce burden on curators.

Talos can be run on any sized dataset, and has been used from single-sample panel analysis to multi-thousand sample
genome callsets. There are benefits to running larger callsets (Talos shows better than linear scaling with increased
sample counts), but analysis is per-family unit, so cohort size doesn't affect results.

## Performance

Talos has been benchmarked on a range of clinical and research datasets, with a publication near completion. In brief,
Talos returns a mean of one variant per proband per family, and 2 variants per singleton, with approximately 25% of all
returned variants, and around 40% of variants flagged for special attention deemed clinically relevant.

In terms of hardware consumption, the demonstration workflow using toy data embedded in this repository takes roughly
2 minutes to run the annotation workflow, and 5 minutes to run Talos. This is all using 2 cpu cores and 8GB memory.

Real world usage will vary wildly depending on the size of the input data, whether merging of VCFs is done within or
prior to the workflow, and whether you're using Exomes or Genomes. I appreciate that this is not very useful, but we
would love to get better usage information to improve this. If you run Talos, please use the `-with-report` command line
parameter, and let us know some details about your run. Sharing this file with us (via email) would be fantastic
information to help us understand and improve the performance of Talos.

### Workflows

Included here are two reference workflow implementations using NextFlow. To run these, NextFlow will have to be installed
on the machine operating the workflow, see [here](https://www.nextflow.io/docs/latest/install.html) for instructions.
The workflow itself runs in a containerised manner using a Docker image built on the [Dockerfile](Dockerfile) at the root
of this repository. This Docker image contains Talos and all its dependencies, plus BCFtools (used for merging and
consequence annotation) and [Echtvar](https://github.com/brentp/echtvar) (used to rapidly apply population frequencies).

- [Annotation](nextflow/annotation.nf): Starting from single or multisample VCFs, this workflow annotates and reformats
  the variant data into a Talos-ready starting point.
- [Talos](nextflow/talos.nf): A full Talos analysis, starting from the annotated data in MatrixTable form, and running
  through to report generation.

> **_NOTE:_**  Note the tag of the dockerfile in this command is kept in sync with the package version and config
> setting. If you apply another tag you'll have to make the corresponding change in the nextflow config files.

```commandline
docker build -t talos:7.1.0 .
```

The [individual Nextflow Modules](nextflow/modules) describe each step of the pipeline, and could be reimplemented in
any other framework. We'd be glad to discuss specific implementations for your use case.

### Example Inputs

A range of stub files have been provided to demonstrate each workflow. This includes a trio of individual VCFs, and a
corresponding pedigree as a toy example. In addition to these stub files you'll need some larger files which are not
economical to store in this repository.
These input files are expected in the `large_files` directory, which in the [default annotation config](nextflow/annotation.config) is at the root of this repository. The corresponding [`README.md`](large_files/README.md) file contains a list of files that need to be downloaded.

All nextflow configuration parameters can be overridden at runtime - to supply an alternate location for the 'large_files' directory use `--large_files XXX` as a CLI parameter when starting the workflow. To alter the output location for the processed annotation files, alter the `--processed_annotations YYY` parameter.

Once these files are downloaded, and a Docker file is built from the root of this repository, you can run the workflows with these commands (optional parameters in square brackets):

```commandline
nextflow -c nextflow/annotation.config run nextflow/annotation.nf [--large_files XXX] [--processed_annotations YYY] [--merged_vcf ZZZ]

and

nextflow -c nextflow/talos.config run nextflow/talos.nf --matrix_tar nextflow/cohort_outputs/cohort.mt.tar
```

For the AlphaMissense, MANE, and Ensembl GFF file, the annotation workflow procesess the raw data into a format expected in the workflow. This only needs to be done once, and subsequent runs of the workflow will detect and reuse the prior outputs, with faster runtimes.

## Real Data

For real data, you'll need to provide your own VCFs, and a pedigree file. Optionally you can supply a Phenopackets file, used to supply additional phenotypic information for participants. This is used to generate a personalised gene panel for each family in the analysis. If this is not provided, an empty phenopacket
object will be generated from the pedigree, and can be updated for subsequent analysis if and when phenotypic information becomes available.

A newly supported argument to the annotation workflow is `--merged_vcf`. This can be used to provide a single merged VCF
file as input, instead of merging individual VCFs. This will be beneficial if you already have merged data. To ensure uniformity the
merged VCF will be passed through a variant normalisation and region filtering step, which can take some time depending on its size.

---

Talos analysis consists of a few main phases:

1. [optional] Using HPO terms in metadata to identify phenotype-matched gene panels of interest
   - If no phenotype terms were provided, analysis will just default to the base (Mendeliome) gene panel
2. Querying PanelApp to find the genes of interest for this analysis run
   - Using any phenotype-matched panels overlaid on a base Mendelian gene panel
3. Variant Filtering & Categorisation
   - Talos applies several filtering criteria to remove common or benign variants, then applies a number of filtering
     categories to select variants which pass a decision tree of criteria.
   - If a variant passes all criteria of a category, it is labelled with that category.
   - If a variant passes multiple categories, it is labelled with all applicable categories.
   - Once all categories have been applied, any un-categorised variants are removed.
4. Mode of Inheritance (MOI) Checking
   - For each remaining variant, we check if the variant's presence in members of the family is consistent with its MOI
     - This includes both individual variants and compound-heterozygotes between multiple different variants.
   - If the variant is consistent with the MOI, the variant is retained for a final report.
5. [optional] HPO Term Matching
   - If phenotypes are provided, we flag any variants which seem well matched to their families to prioritise analysis
6. Report Generation
   - An HTML report is generated for the cohort as a whole, and separately for each family, detailing the variants which
     passed all filtering criteria

### Input Data

To run Talos you will need:

1. Variant data, annotated and formatted by the enclosed annotation process. To use alternative entrypoints (e.g. an existing MatrixTable or annotated VCF, please get in touch - we can help with adapting the data, but it's unlikely to be in the expected format)
2. ClinVar data as generated by ClinvArbitration. This is available from the [ClinvArbitration Release Page](https://github.com/populationgenomics/ClinvArbitration/releases), or can be generated using the code and process described in the [ClinvArbitration repository](https://github.com/populationgenomics/ClinvArbitration).
3. A pedigree file, describing the pedigree of the participants in the study ([Pedigree Reference](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format))
4. A TOML file containing the configuration for the analysis. [example_config.toml](src/talos/example_config.toml) is a good starting point, with comments explaining each modifiable parameter. Changes you may wish to make to tailor analysis to the cohort under test are using `forced_panels` to involve additional gene panels
   in the analysis, removing or extending the `require_pheno_match` list, which would mask noisy variants from the base panel, and `forbidden_genes` to remove genes from the analysis completely.
6. [Optional] A JSON file detailing a GA4GH compliant Cohort ([see reference](https://phenopacket-schema.readthedocs.io/en/latest/cohort.html)). This should contain Phenotypic Features for all relevant participants using HPO terms. If provided this improves matching of panels to participants. If you previously generated the
   'extended PED' file in place of the phenopackets file, there is a conversion script [here: convert_ePED_to_phenopackets.py](src/talos/ConvertPedToPhenopackets.py) which will convert the extended PED file to a Phenopackets file and regular pedigree. If this is not provided at all, an empty phenopacket file will be derived
   from the pedigree.

## Usage

Talos consists of the following components:

- `MakePhenopackets` - This is a CPG-specific implementation for generating a Cohort/Phenopacket file. It can serve as a template for generating a compliant Phenopackets input file.
- `DownloadPanelApp` - intended to run once per month, this process downloads a simplified representation of the whole of panelapp at a given timepoint. This can be run outside of the main workflow if you are deploying to an offline compute environment.
- `UnifiedPanelAppParser` - Takes the Downloaded panelapp data and the phenopacket file, and generates the region of interest for the analysis. If a HPO ontology file is provided (see [large_files](large_files/README.md)) it will run phenotypic matching to tailor the region of interest to each proband.
- `RunHailFiltering` - Takes the MatrixTable of Variants, the Pedigree file, the panel data from `UnifiedPanelAppParser`, and both ClinVar tables, filters the variants in the MatrixTable, and labels them with categories of interest.
  This is the most resource-intensive step of the pipeline, but even on 400+GB datasets it has been run successfully on a 8-core, 16GB RAM VM.
- `RunHailFilteringSv` - Takes a MatrixTable of Structural Variants, the Pedigree file, the panel data from `UnifiedPanelAppParser`, filters the variants in the MatrixTable, and labels them with categories of interest.
- `ValidateMOI` - Takes the result of `RunHailFiltering`, optionally an SV result from `RunHailFilteringSV`, the Pedigree, and panel data from `UnifiedPanelAppParser`. Checks each categorised variant to determine whether the MOI associated with the relevant gene fits within the family structure where it occurs.
- `HPOFlagging` - Takes the results of `ValidateMOI`, and uses semsimian to test whether the HPO term(s) associated with the gene matches the HPO term(s) associated with the participant. Used in escalated variant prioritisation by flagging where a phenotypic match has been detected.
- `CreateTalosHTML` - Generates a report from the results of the `ValidateMOI`.
- `MinimiseOutputForSeqr` - Parses the result of `ValidateMOI`, generates a file for ingestion by Seqr.

## Categories

![CatDiagram](design_docs/images/Categories.png)

This is a highly simplified representation of the categories currently implemented.

## ClinvArbitration

See the companion [ClinvArbitration](https://github.com/populationgenomics/ClinvArbitration) repository for more details.

Talos uses ClinVar submissions to determine if a variant has been previously reported as Pathogenic or Likely Pathogenic. For the purpose of a Talos analysis, when multiple submissions for the same variant have conflicting classifications, we would prefer to favour the classifications provided by high-quality recent submitters. To enable this, we have ClinvArbitration, a re-aggregation of ClinVar submissions uses altered heuristics favouring decisions. This avoids the default logic used by clinvar which cautiously assigns "conflicting interpretations of pathogenicity" unless there is perfect harmony among all submissions.

In addition to providing the top-line rating, ClinvArbitration also re-indexes ClinVar variants based on the transcript and codon they alter. Talos uses this table to identify any missense variants impacting a codon previously reported to be altered by a pathogenic variant in ClinVar and assign the `PM5` evidence category.

ClinvArbitration is used to re-process ClinVar releases periodically and the ready-to-use results are available via the ClinvArbitration repository's Release page.

## Exomiser Integration

[Exomiser](https://github.com/exomiser/Exomiser) is a similar tool in this space, running variant prioritisation based on variant annotations and matching genes to participant phenotypes. In our benchmarking cohorts we've found that performance of Talos and Exomiser are comparable when used independently, but when the two analyses are combined, they perform better than either in isolation. To that end we have facilitated the integration of Exomiser variants into Talos analyses. If you have previously generated Exomiser (v14) results, the script [CombineExomiserVariantTsvs.py](src/talos/AggregateExomiserVariantTsvs.py) will aggregate the per-proband TSVs into a single JSON file, and write the results as a Hail Table. This Hail Table can be passed to `RunHailFiltering`, with variants prioritised by Exomiser annotated in as another data source. In the HTML report, an `Exomiser` column will display the rank and Mode of Inheritance (MOI) of the variant in the Exomiser reports.

This preparatory script expects that the Exomiser Variant analysis TSVs have the naming convention `PROBAND.variant.tsv`, where PROBAND is the name of the proband in the cohort/VCF/MatrixTable.

## Reanalysis

The heart of Talos' utility is in re-analysis, by bootstrapping from previous analyses. Where possible each run consults the history from the previous run, determining whether each variant has been seen before, and if so, whether evidence has evolved. Each run adds the incremental content, and re-saves the history.

The final report contains a `first_seen` date for each variant, along with an `evidence_last_updated` date which indicates the most recent date that the evidence changed (new category labels were applied). By filtering on either of these dates, analysts can view only the incremental variants new in each round, or variants where the evidence has changed.
