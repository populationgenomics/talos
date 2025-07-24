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
include { FilterVcfToBedWithBcftools } from './modules/annotation/FilterVcfToBedWithBcftools/main'
include { MakeSitesOnlyVcfWithBcftools } from './modules/annotation/MakeSitesOnlyVcfWithBcftools/main'
include { MergeVcfsWithBcftools } from './modules/annotation/MergeVcfsWithBcftools/main'
include { ParseAlphaMissenseIntoHt } from './modules/annotation/ParseAlphaMissenseIntoHt/main'
include { ParseManeIntoJson } from './modules/annotation/ParseManeIntoJson/main'
include { ReformatAnnotatedVcfIntoHailTable } from './modules/annotation/ReformatAnnotatedVcfIntoHailTable/main'
include { TransferAnnotationsToMatrixTable } from './modules/annotation/TransferAnnotationsToMatrixTable/main'

workflow {

    // populate ref genome input channel
    ch_ref_genome = channel.fromPath(
    	params.ref_genome,
    	checkIfExists: true,
	)

    // generate the AlphaMissense HT - long running, stored in a separate folder
    // read in as a channel if this was already generated
    if (file(params.alphamissense_tar).exists()) {
        ch_alphamissense_table = channel.fromPath(
        	params.alphamissense_tar,
        	checkIfExists: true
		)
    }
    else {
    	ch_alphamissense_tsv = channel.fromPath(
    		params.alphamissense_tsv,
    		checkIfExists: true
		)
        ParseAlphaMissenseIntoHt(ch_alphamissense_tsv)
        ch_alphamissense_table = ParseAlphaMissenseIntoHt.out
    }

    // generate the Region-of-interest BED file from Ensembl GFF3
    // generates a per-gene BED file with ID annotations
    // and a overlap-merged version of the same for more efficient region filtering
    ch_gff = channel.fromPath(params.ensembl_gff, checkIfExists: true)
    if (file(params.ensembl_bed).exists() && file(params.ensembl_merged_bed).exists()) {
    	ch_bed = channel.fromPath(
    		params.ensembl_bed,
    		checkIfExists: true,
		)
    	ch_merged_bed = channel.fromPath(
    		params.ensembl_merged_bed,
    		checkIfExists: true,
		)
    }
    else {
    	CreateRoiFromGff3(ch_gff)
    	ch_bed = CreateRoiFromGff3.out.bed
    	ch_merged_bed = CreateRoiFromGff3.out.merged_bed
	}

	// if a merged VCF is provided, don't implement a manual merge - start from an externally completed dataset
	if (file(params.merged_vcf).exists()) {
		ch_merged_vcf = channel.fromPath(params.merged_vcf, checkIfExists: true)
		ch_merged_index = channel.fromPath("${params.merged_vcf}.tbi", checkIfExists: true)

		FilterVcfToBedWithBcftools(
			ch_merged_vcf,
			ch_merged_index,
			ch_merged_bed,
			ch_ref_genome,
		)
		ch_merged_tuple = FilterVcfToBedWithBcftools.out
	}
	else {
		ch_vcfs = channel.fromPath("${params.input_vcf_dir}/*.vcf.gz")
		ch_tbis = channel.fromPath("${params.input_vcf_dir}/*.vcf.gz").map{ it -> file("${it}.tbi") }
		MergeVcfsWithBcftools(
			ch_vcfs.collect(),
			ch_tbis.collect(),
			ch_merged_bed,
			ch_ref_genome,
		)
		ch_merged_tuple = MergeVcfsWithBcftools.out
	}

    // create a sites-only version of this VCF, just to pass less data around when annotating
    MakeSitesOnlyVcfWithBcftools(
        ch_merged_tuple
    )

    // read the whole-genome Zip file as an input channel
    ch_gnomad_zip = channel.fromPath(
		params.gnomad_zip,
		checkIfExists: true
    )

    // and apply the annotation using echtvar
    AnnotateGnomadAfWithEchtvar(
        MakeSitesOnlyVcfWithBcftools.out,
        ch_gnomad_zip,
    )

    // annotate transcript consequences with bcftools csq
    AnnotateCsqWithBcftools(
        AnnotateGnomadAfWithEchtvar.out,
        ch_gff,
        ch_ref_genome,
    )

    // pull and parse the MANE data into a Hail Table
    if (file(params.mane_json).exists()) {
    	ch_mane = channel.fromPath(
			params.mane_json,
			checkIfExists: true
		)
    }
    else {
    	ch_mane_summary = channel.fromPath(
    		params.mane,
    		checkIfExists: true
		)
    	ParseManeIntoJson(ch_mane_summary)
    	ch_mane = ParseManeIntoJson.out.json
    }

    // reformat the annotations in the VCF, retain as a Hail Table
    ReformatAnnotatedVcfIntoHailTable(
        AnnotateCsqWithBcftools.out,
        ch_alphamissense_table,
        ch_bed,
        ch_mane,
    )

    // combine the join-VCF and annotations as a HailTable
    TransferAnnotationsToMatrixTable(
        ReformatAnnotatedVcfIntoHailTable.out,
        ch_merged_tuple,
    )
}
