#!/usr/bin/env nextflow

/*
This workflow is the annotation process for the Talos pipeline.

It requires a single VCF, which can be single- or multi-sample. This is then annotated and reformatted for Talos.
The specific annotations are:

- gnomAD v4.1 frequencies and alphamissense annotations, applied to the joint VCF using echtvar
- Transcript consequences, using BCFtools annotate
- MANE trancript IDs and corresponding ENSP IDs, applied using Hail

TOOD: Planned extension - multiple VCFs can be provided, instead of using the VCF splitting workflow.
*/

nextflow.enable.dsl=2

include { AnnotateCsqWithBcftools } from './modules/annotation/AnnotateCsqWithBcftools/main'
include { AnnotatedVcfIntoMatrixTable } from './modules/annotation/AnnotatedVcfIntoMatrixTable/main'
include { AnnotateWithEchtvar } from './modules/annotation/AnnotateWithEchtvar/main'
include { NormaliseAndRegionFilterVcf } from './modules/annotation/NormaliseAndRegionFilterVcf/main'
include { SplitVcf } from './modules/annotation/SplitVcf/main'


workflow {
    main:
    // populate various input channels - these are downloaded by the large_files/gather_file.sh script
    ch_gff = channel.fromPath(params.ensembl_gff, checkIfExists: true)
    ch_gnomad_zip = channel.fromPath(params.gnomad_zip, checkIfExists: true)
    ch_ref_genome = channel.fromPath(params.ref_genome, checkIfExists: true)

    // check that the JSON-format MANE data exists already, or prompt to run the prep workflow
    if (!file(params.mane_json).exists()) {
        println "MANE JSON not available, please run the Talos Prep workflow (talos_preparation.nf)"
        exit 1
    }

    ch_mane = channel.fromPath(params.mane_json, checkIfExists: true)

    // Read the AlphaMissense HT as a channel, or prompt for generation using the prep workflow
    if (!file(params.alphamissense_zip).exists()) {
        println "AlphaMissense data must be encoded for echtvar, run the Talos Prep workflow (talos_preparation.nf)"
        exit 1
    }

    ch_alphamissense_zip = channel.fromPath(params.alphamissense_zip, checkIfExists: true)

    // check the ensembl BED file has been generated
    if (!file(params.ensembl_bed).exists()) {
        println "Region-Of-Interest BED file has not been prepared, run the Talos Prep workflow (talos_preparation.nf)"
        exit 1
    }
    ch_bed = channel.fromPath(params.ensembl_bed, checkIfExists: true)
    ch_merged_bed = channel.fromPath(params.ensembl_merged_bed, checkIfExists: true)

    // see if sharded VCFs were provided
    if (params.shards != null) {
        ch_vcfs = Channel.fromPath("${params.shards}/*.${params.input_vcf_extension}")
    } else {
        ch_vcf = channel.fromPath(params.vcf, checkIfExists: true)
        // decide whether to split and parallelise, or run as a single operation
        // if config value is absent completely, skip this step
        if ((params.vcf_split_n ?: 0) > 0) {
            SplitVcf(ch_vcf)
            ch_vcfs = SplitVcf.out.flatten()
        } else {
            ch_vcfs = ch_vcf
        }
    }

	NormaliseAndRegionFilterVcf(
        ch_vcfs,
        ch_merged_bed.first(),
        ch_ref_genome.first(),
    )

	AnnotateWithEchtvar(
        NormaliseAndRegionFilterVcf.out,
        ch_gnomad_zip.first(),
        ch_alphamissense_zip.first(),
    )

    // annotate transcript consequences with bcftools csq
    AnnotateCsqWithBcftools(
        AnnotateWithEchtvar.out,
        ch_gff.first(),
        ch_ref_genome.first(),
    )

    // reformat the annotations in the VCF, generate a Hail MatrixTable
    AnnotatedVcfIntoMatrixTable(
        AnnotateCsqWithBcftools.out,
        ch_bed.first(),
        ch_mane.first(),
    )

}
