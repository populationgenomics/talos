# SV Annotation Module

A Nextflow workflow which annotates a joint-called Structural Variant VCF with:

- **gene consequences**, using [GATK `SVAnnotate`](https://gatk.broadinstitute.org/hc/en-us/articles/13832752531611-SVAnnotate) against the MANE GTF
- **population allele frequencies**, using [SVAFotate](https://github.com/fakedrtom/SVAFotate) against gnomAD-SV

This is a re-implementation of the same core steps we (CPG) use internally. Our internal usage centres around the
GATK-SV workflow, and our CPG-Flow wrapped implementation of it. The terminal stage of this workflow is [Annotation](https://github.com/populationgenomics/cpg-flow-gatk-sv/blob/main/src/cpg_flow_gatk_sv/multisample_workflow.py#L701),
which is done using GATK's SvAnnotate tool for consequence prediction, and a complex interval-overlap-resolution step
to match variants to gnomAD frequencies.

Instead of re-implementing the exact process here, I've split the annotation into two phases:

  - Consequence: handled using SVAnnotate, and exact replica of the GATK-SV process
  - Pop.Freq: handled using [SVAFotate](https://github.com/fakedrtom/SVAFotate)

These two steps, and pre-processing of relevant input files, run for each cohort whose input TSV row has an `sv` VCF,
with the same core conceit as small variants and Mito data - a single joint-called VCF should contain the whole
group of Samples being processed, which should also match the Pedigree and Small-variant data.

A separate sub-workflow, `SV_ANNOTATION` has been created to handle these steps. `SV_ANNOTATION` publishes an annotated
VCF per cohort, `RunSmallFilteringSv` filters and labels it with `CategoryBooleanSV1`, and `ValidateMOI` folds the result
into the report.

To utilise this functionality, add a `sv` column to the input TSV, pointing to a SV joint-call. This has been tested on
the output of GATK-SV (multiple variant callers, with resolved calls) and GATK's gCNV. The decision is made per row:
cohorts with an empty or absent `sv` cell skip SV annotation and filtering.

## Pre-annotated SV VCFs

If a cohort's SV VCF has already been annotated (e.g. the output of GATK-SV's own annotation stage), put it in the
`sv_annotated` column instead of `sv`. It skips `SV_ANNOTATION` entirely and goes straight to `RunSvFiltering`.
The VCF must:

- be bgzipped, with a `.tbi` index alongside it (`<vcf>.tbi`)
- carry SVAnnotate's consequence fields (`PREDICTED_LOF` etc.)
- carry population frequencies as `<prefix>_sv_AF` and `<prefix>_sv_SVID`, where `<prefix>` is the cohort config's
  `RunSvFiltering.gnomad_population`. These fields are not renamed, so set the config to match the VCF

A row can set `sv` or `sv_annotated`, not both; the run fails before any task starts if one does.

The SV reference data (`svafotate_bed`, `mane_gtf`, `svannotate_noncoding_bed`, all downloaded by
`large_files/gather_files.sh`) is required for every run, whether or not any cohort has SV data. `ref_dict` is only
generated when there is SV annotation to do.
