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

include { AnnotateCsqWithBcftools } from './modules/annotation/AnnotateCsqWithBcftools/main'
include { AnnotateGnomadAfWithEchtvar } from './modules/annotation/AnnotateGnomadAfWithEchtvar/main'
include { CreateRoiFromGff3 } from './modules/annotation/CreateRoiFromGff3/main'
include { LocaliseAlphamissenseWithWget } from './modules/annotation/LocaliseAlphamissenseWithWget/main'
include { MakeSitesOnlyVcfWithBcftools } from './modules/annotation/MakeSitesOnlyVcfWithBcftools/main'
include { MergeVcfsWithBcftools } from './modules/annotation/MergeVcfsWithBcftools/main'
include { ParseAlphaMissenseIntoHt } from './modules/annotation/ParseAlphaMissenseIntoHt/main'
include { ParseManeIntoJson } from './modules/annotation/ParseManeIntoJson/main'
include { ReformatAnnotatedVcfIntoHailTable } from './modules/annotation/ReformatAnnotatedVcfIntoHailTable/main'
include { TransferAnnotationsToMatrixTable } from './modules/annotation/TransferAnnotationsToMatrixTable/main'

workflow {

    // populate input channels - VCFs, reference genome
    ch_vcfs = channel.fromPath(params.input_vcfs)
    ch_tbis = channel.fromPath(params.input_vcfs).map{ it -> file("${it}.tbi") }
    ch_ref_genome = channel.fromPath(params.ref_genome)

    // generate the AlphaMissense HT - long running, stored in a separate folder
    // read in as a channel if this was already generated
    if (file(params.alphamissense_output).exists()) {
        ch_alphamissense_table = channel.fromPath(params.alphamissense_output)
    }
    else {
        LocaliseAlphamissenseWithWget()
        ParseAlphaMissenseIntoHt(LocaliseAlphamissenseWithWget.out)
        ch_alphamissense_table = ParseAlphaMissenseIntoHt.out
    }

    // generate the Region-of-interest BED file from Ensembl GFF3
    // generates a per-gene BED file with ID annotations
    // and a overlap-merged version of the same for more efficient region filtering
    CreateRoiFromGff3()

    MergeVcfsWithBcftools(
        ch_vcfs.collect(),
        ch_tbis.collect(),
        CreateRoiFromGff3.out.merged_bed,
    )

    // create a sites-only version of this VCF, just to pass less data around when annotating
    MakeSitesOnlyVcfWithBcftools(
        MergeVcfsWithBcftools.out
    )

    // read the whole-genome Zip file as an input channel
    ch_gnomad_zip = channel.fromPath(params.gnomad_zip)

    // and apply the annotation using echtvar
    AnnotateGnomadAfWithEchtvar(
        MakeSitesOnlyVcfWithBcftools.out,
        ch_gnomad_zip,
    )

    // annotate transcript consequences with bcftools csq
    AnnotateCsqWithBcftools(
        AnnotateGnomadAfWithEchtvar.out,
        CreateRoiFromGff3.out.gff3,
        ch_ref_genome
    )

    // pull and parse the MANE data into a Hail Table
    ParseManeIntoJson()

    // reformat the annotations in the VCF, retain as a Hail Table
    ReformatAnnotatedVcfIntoHailTable(
        AnnotateCsqWithBcftools.out,
        ch_alphamissense_table,
        CreateRoiFromGff3.out.bed,
        ParseManeIntoJson.out.json
    )

    // combine the join-VCF and annotations as a HailTable
    TransferAnnotationsToMatrixTable(
        ReformatAnnotatedVcfIntoHailTable.out,
        MergeVcfsWithBcftools.out,
    )
}
