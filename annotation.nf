#!/usr/bin/env nextflow

/*
This workflow is a minimised annotation process for the Talos pipeline.

It takes multiple VCFs, merges them into a single VCF, and annotates the merged VCF with relevant annotations.

The specific annotations are:
- gnomAD v4.1 frequencies, applied to the joint VCF using echtvar
- Transcript consequences, using BCFtools annotate
- AlphaMissense, applied using Hail
- MANE trancript IDs and corresponding proteins, applied using Hail
*/

nextflow.enable.dsl=2

include { alphamissense_to_ht } from './modules/annotation/alphamissense_to_ht/main'
include { csq_annotation } from './modules/annotation/csq_annotation/main'
include { generate_roi } from './modules/annotation/generate_roi/main'
include { gnomad_annotate } from './modules/annotation/gnomad_annotate/main'
include { localise_alphamissense } from './modules/annotation/localise_alphamissense/main'
include { merge_vcfs } from './modules/annotation/merge_vcfs/main'
include { pull_and_parse_mane } from './modules/annotation/pull_and_parse_mane/main'
include { vcf_to_mt } from './modules/annotation/vcf_to_mt/main'

workflow {
    // generate the AlphaMissense HT - long running, stored in a separate folder
    if (file(params.alphamissense_output).exists()) {
        alphamissense_table = channel.fromPath(params.alphamissense_output)
    }
    else {
        localise_alphamissense()
        alphamissense_to_ht(localise_alphamissense.out)
        alphamissense_table = alphamissense_to_ht.out
    }

    // generate the gene region file
    generate_roi()

    // get all the VCFs
    vcfs = channel.fromPath(params.input_vcfs)
    tbis = channel.fromPath(params.input_vcfs).map{ it -> file("${it}.tbi") }

    merge_vcfs(
        vcfs.collect(),
        tbis.collect(),
        generate_roi.out.bed,
    )

    // now get the annotation channels (one zip per chromosome)
    annotation_zip = channel.fromPath(params.gnomad_zip)

    // and apply the annotation using echtvar
    gnomad_annotate(merge_vcfs.out, annotation_zip)

    // bcftools csq
    ref_genome = channel.fromPath(params.ref_genome)
    csq_annotation(gnomad_annotate.out, generate_roi.out.gff3, ref_genome)

    // pull and parse the MANE data into a Hail Table
    pull_and_parse_mane()

    // now what
    vcf_to_mt(csq_annotation.out, alphamissense_table, generate_roi.out.bed, pull_and_parse_mane.out)
}


