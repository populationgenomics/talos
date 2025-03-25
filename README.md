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

### Example Workflows

Included here are two reference workflow implementation using NextFlow:

- [Annotation](nextflow/annotation.nf): Starting from single or multisample VCFs, this workflow annotates and reformats
    the variant data into a Talos-ready starting point.
- [Talos](nextflow/talos.nf): A full Talos analysis, starting from the annotated data in MatrixTable form, and running
    through to report generation.

These workflows can be executed in full using a Docker image built on the [Dockerfile](Dockerfile) at the root of this
repository. This Docker image contains Talos and all its dependencies, plus BCFtools (used for merging and consequence
annotation) and [Echtvar](https://github.com/brentp/echtvar) (used to rapidly apply population frequencies). 

> **_NOTE:_**  Note the tag of the dockerfile in this command is kept in sync with the package version and config 
> setting. If you apply another tag you'll have to make the corresponding change in the nextflow config files.

```commandline
docker buildx build -t talos:7.0.1 .
```

> **_NOTE:_**  Note the tag of the dockerfile in this command is kept in sync with the package version and config
> setting. If you apply another tag you'll have to make the corresponding change in the nextflow config files.

```commandline
docker buildx build -t talos:7.0.1 .
```

The [individual Nextflow Modules](nextflow/modules) describe each step of the pipeline, and could be reimplemented in
any other framework. We'd be glad to discuss specific implementations for your use case.

### Input Data

A range of stub files have been provided to demonstrate each workflow. This includes a trio of individual VCFs, and a
corresponding pedigree as a toy example. In addition to these stub files you'll need some larger files which are not
economical to store in this repository. For the annotation workflow:

1. A reference genome matching your input data, in FASTA format
2. An echtvar reference file - we have a pre-generated file we are able to share, please get in touch

For the Talos workflow:

1. `genes_to_phenotypes.txt` from https://hpo.jax.org/data/annotations
2. `hp.obo` from https://hpo.jax.org/data/ontology
3. `phenio.db.gz` from https://data.monarchinitiative.org/monarch-kg/latest/phenio.db.gz

These are expected in the [large_files](large_files) directory at the root of this repository.

Once these files are downloaded, and a Docker file is built from the root of this repository, you can run the workflows
with:

```commandline
nextflow -c nextflow/annotation.config run nextflow/annotation.nf

and

nextflow -c nextflow/talos.config run nextflow/talos.nf --matrix_tar nextflow/cohort_outputs/cohort.mt.tar
```

The first time the annotation workflow runs will be slower - large AlphaMissense data is downloaded from zenodo, but
this only needs to be done on the first run, and a couple of ancillary files are generated to use in annotation.

### Real Data

For real data, you'll need to provide your own VCFs, and a pedigree file. Optionally you can supply a Phenopackets file,
used to supply additional phenotypic information for participants. This is used to generate a personalised gene panel
for each family in the analysis.

---

Talos analysis consists of a few main phases:

1. [optional] Using HPO terms in metadata to identify phenotype-matched gene panels of interest
    * If no phenotype terms were provided, analysis will just default to the base (Mendeliome) gene panel
2. Querying PanelApp to find the genes of interest for this analysis run
   * Using any phenotype-matched panels overlaid on a base Mendelian gene panel
3. Variant Filtering & Categorisation
   * Talos applies several filtering criteria to remove common or benign variants, then applies a number of filtering
        categories to select variants which pass a decision tree of criteria.
   * If a variant passes all criteria of a category, it is labelled with that category.
   * If a variant passes multiple categories, it is labelled with all applicable categories.
   * Once all categories have been applied, any un-categorised variants are removed.
4. Mode of Inheritance (MOI) Checking
   * For each remaining variant, we check if the variant's presence in members of the family is consistent with its MOI
     * This includes both individual variants and compound-heterozygotes between multiple different variants.
   * If the variant is consistent with the MOI, the variant is retained for a final report.
5. [optional] HPO Term Matching
   * If phenotypes are provided, we flag any variants which seem well matched to their families to prioritise analysis
6. Report Generation
   * An HTML report is generated for the cohort as a whole, and separately for each family, detailing the variants which
       passed all filtering criteria


## Input Data

To run Talos you will need:

1. Variant data, annotated with VEP. The input can be provided as a Hail MatrixTable or as a multisample VCF

    * Talos uses Hail Query, a PySpark-based query engine, to perform highly parallelised analysis. This requires variants to be stored using the Hail MatrixTable format. If your current workflow uses hail, a MatrixTable can be provided directly as an input.
    * Alternatively a VEP-annotated multi-sample VCF can be provided as input. An additional pre-processing step will convert the VCF to a MatrixTable at run time.
    * Talos is intended to run once per-cohort, not once per cohort. Variant calls from all families/individuals in a cohort should be merged into a single multi-sample file prior to processing with Talos.

2. ClinVar data as generated by ClinvArbitration, both the `clinvar_decisions` and `clinvar_pm5` Hail Tables. This is
   available from the [ClinvArbitration Release Page](https://github.com/populationgenomics/ClinvArbitration/releases), or can be generated using the code and process described in the [ClinvArbitration repository](https://github.com/populationgenomics/ClinvArbitration).
3. A pedigree file, describing the pedigree of the participants in the study ([Pedigree Reference](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format))
4. A TOML file containing the configuration for the analysis. [example_config.toml](src/talos/example_config.toml) is a
   good starting point, with comments explaining each modifiable parameter. Changes you may wish to make to tailor
   analysis to the cohort under test are using `forced_panels` to involve additional gene panels in the analysis,
   removing or extending the `require_pheno_match` list, which would mask noisy variants from the base panel, and
   `forbidden_genes` to remove genes from the analysis completely.
5. [Optional] A JSON file detailing a GA4GH compliant Cohort ([see reference](https://phenopacket-schema.readthedocs.io/en/latest/cohort.html)). This should contain Phenotypic Features for all relevant participants using HPO terms. If provided this improves matching of panels to participants. If you previously generated the 'extended PED' file in place of the phenopackets file, there is a conversion script [here: convert_ePED_to_phenopackets.py](src/talos/ConvertPedToPhenopackets.py) which will convert the extended PED file to a Phenopackets file and regular pedigree.

## Usage

Talos consists of the following components:

- `MakePhenopackets` - This is a CPG-specific implementation for generating a Cohort/Phenopacket file. It can serve as a template for generating a compliant Phenopackets input file.
- `GeneratePanelData` - [optional] Phenopacket file, and generates a per-participant list of panels to be used for this analysis, writing the result as a JSON. This also requires a local copy of the HPO ontology, downloadable from [here](http://purl.obolibrary.org/obo/hp.obo).
- `QueryPanelapp` - Takes the output of `GeneratePanelData`, or None if no PhenoPacket file was provided, Queries PanelApp for the panels selected for the cohort, and writes the result as a JSON.
- `FindGeneSymbolMap` - Uses the output of `QueryPanelapp` to find the gene symbol for each gene ID via Ensembl's REST API.
- `RunHailFiltering` - Takes the MatrixTable of Variants, the Pedigree file, the panel data from `QueryPanelapp`, and
  both ClinVar tables, filters the variants in the MatrixTable, and labels them with categories of interest. This is the most resource-intensive step of the pipeline, but even on 400+GB datasets it has been run successfully on a 8-core, 16GB RAM VM.
- `RunHailFilteringSv` - Takes a MatrixTable of Structural Variants, the Pedigree file, the panel data
  from `QueryPanelapp`, filters the variants in the MatrixTable, and labels them with categories of interest.
- `ValidateMOI` - Takes the result of `RunHailFiltering`, optionally one or more SV result from `RunHailFilteringSV`,
  the Pedigree, and panel data from `QueryPanelapp`. Checks each categorised variant to determine whether the MOI
  associated with the relevant gene fits within the family structure where it occurs. Generates a JSON file from all
  variants which pass the MOI tests.
- `HPOFlagging` - Takes the results of `ValidateMOI`, and uses semsimian to test whether the HPO term(s) associated with the gene matches the HPO term(s) associated with the participant.
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
