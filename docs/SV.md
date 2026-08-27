# 🏗️ SV Annotation Module

Since version 11.1.x Talos has been capable of processing SV data. This is done using a sub-workflow, [sv_annotation](nextflow/sv_annotation.nf), which first annotates with consequences and population frequncies, then runs a filter & labelling operation. Selected SVs are considered alongside the small variant data, i.e. a deletion and a stop-gained in a biallelic gene could be flagged for further investigation.

A Nextflow workflow which annotates a joint-called Structural Variant VCF with:

- **gene consequences**, using [GATK `SVAnnotate`](https://gatk.broadinstitute.org/hc/en-us/articles/13832752531611-SVAnnotate) against the MANE GTF
- **population allele frequencies**, using [SVAFotate](https://github.com/fakedrtom/SVAFotate) against gnomAD-SV

This is a re-implementation of the same core steps we (CPG) use internally. Our internal usage centres around the
GATK-SV workflow, and our CPG-Flow wrapped implementation of it. The terminal stage of this workflow is [Annotation](https://github.com/populationgenomics/cpg-flow-gatk-sv/blob/main/src/cpg_flow_gatk_sv/multisample_workflow.py#L701),
which is done using GATK's SvAnnotate tool for consequence prediction, and a complex interval-overlap-resolution step
to match variants to gnomAD frequencies.

Instead of re-implementing the exact process here, I've split the annotation into two phases:

  - Consequence: handled using SVAnnotate, an exact replica of the GATK-SV process
  - Pop.Freq: handled using [SVAFotate](https://github.com/fakedrtom/SVAFotate)

These two steps, and pre-processing of relevant input files, are engaged only if an SV file is included in the input
TSV file, with the same core conceit as small variants and Mito data - a single joint-called VCF should contain the whole
group of Samples being processed, which should also match the Pedigree and Small-variant data.

A separate sub-workflow, `SV_ANNOTATION` has been created to handle these steps. `SV_ANNOTATION` publishes an annotated
VCF per cohort, `RunHailFilteringSv` filters and labels it with `CategoryBooleanSV1`, and `ValidateMOI` folds the result
into the report.

## Inputs

To provide SV data, generate a VCF containing structural variant calls for the samples in the corresponding small-variant VCF & pedigree. This process has been tested with outputs from GATK's gCNV, and the GATK-SV multi-tool SV calling workflow, though it may barf on some other variant callers/VCF formats - please raise an issue if this affects you. 

This VCF should be added to the input TSV file as a `sv` column. The sample sets do not have to overlap completely, i.e. SV results for only one member of a trio, or half of a cohort is perfectly fine. The samples in the SV and small variant VCFs will be compared against the same pedigree file, and any sample IDs not in that pedigree will be ignored.  

Small-variants remain a core component of Talos, and it is not designed to run exclusively on SV data.

To use the SV workflow you will need to:

1. run the latest version of [large_files/gather_files.sh](large_files/gather_files.sh) to download the SVAFotate resources
2. build the svafotate docker file (`docker build -f docker/SVAFotate_Dockerfile -t svafotate:0.1.0 .`)
3. add your SV data under the `sv` input column in the TSV file

The workflow will detect the presence of the SV input, and trigger the SV annotation workflow.

> **Note:** If the SV column is empty but present, you will need the SVAFotate file to be downloaded. If the `sv` column is omitted completely, Talos will run without checking the existence of any SV-specific reference files.

## Demonstration

This repository contains test data to demonstrate the SV functionality. If you run the test workflow using the input file `nextflow/inputs/test_sv.tsv`, this will trigger the SV sub-workflow, and will generate a report containing both small-variant and SV data.

```bash
nextflow \
  -c nextflow.config \
  run main.nf \
  --input_tsv nextflow/inputs/test_sv.tsv \
  -output-dir <path_to_output_dir>
```
