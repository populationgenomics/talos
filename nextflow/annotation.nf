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
include { AnnotateWithEchtvar } from './modules/annotation/AnnotateWithEchtvar/main'
include { CreateRoiFromGff3 } from './modules/annotation/CreateRoiFromGff3/main'
include { EncodeAlphaMissense } from './modules/annotation/EncodeAlphaMissense/main'
include { FilterVcfToBedWithBcftools } from './modules/annotation/FilterVcfToBedWithBcftools/main'
include { MakeSitesOnlyVcfsFromMt } from './modules/annotation/MakeSitesOnlyVcfsFromMt/main'
include { MakeSitesOnlyVcfWithBcftools } from './modules/annotation/MakeSitesOnlyVcfWithBcftools/main'
include { MergeVcfsWithBcftools } from './modules/annotation/MergeVcfsWithBcftools/main'
include { ParseAlphaMissense } from './modules/annotation/ParseAlphaMissense/main'
include { ParseManeIntoJson } from './modules/annotation/ParseManeIntoJson/main'
include { ReformatAnnotatedVcfIntoHailTable } from './modules/annotation/ReformatAnnotatedVcfIntoHailTable/main'
include { TransferAnnotationsToMatrixTable } from './modules/annotation/TransferAnnotationsToMatrixTable/main'
include { VcfToMatrixTable } from './modules/annotation/VcfToMatrixTable/main'

workflow {

    // populate ref genome input channel
    ch_ref_genome = channel.fromPath(
    	params.ref_genome,
    	checkIfExists: true,
	)

    // generate the AlphaMissense HT - long running, stored in a separate folder
    // read in as a channel if this was already generated
    if (file(params.alphamissense_zip).exists()) {
        ch_alphamissense_zip = channel.fromPath(
        	params.alphamissense_zip,
        	checkIfExists: true
		)
    }
    else {
    	ch_alphamissense_tsv = channel.fromPath(
    		params.alphamissense_tsv,
    		checkIfExists: true
		)
        ParseAlphaMissense(ch_alphamissense_tsv)
        EncodeAlphaMissense(ParseAlphaMissense.out)
        ch_alphamissense_zip = EncodeAlphaMissense.out
    }

    // read the whole-genome Zip file as an input channel
    ch_gnomad_zip = channel.fromPath(
		params.gnomad_zip,
		checkIfExists: true
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

	// check if we're using a MT as input
	if (file("${params.matrix_table}/_SUCCESS").exists()) {
	    ch_mt = channel.fromPath(params.matrix_table, checkIfExists: true)
	}
	else {
	    // if a merged VCF is provided, don't implement a manual merge - start from an externally completed dataset
        if (file(params.merged_vcf).exists()) {
            VcfToMatrixTable(
                channel.fromPath(params.merged_vcf, checkIfExists: true),
                ch_merged_bed,
            )
            ch_mt = VcfToMatrixTable.out
        }
        else {
            ch_vcfs = channel.fromPath("${params.input_vcf_dir}/*.${params.input_vcf_extension}", checkIfExists: true )
            ch_tbis = ch_vcfs.map{ it -> file("${it}.tbi") }
            MergeVcfsWithBcftools(
                ch_vcfs.collect(),
                ch_tbis.collect(),
                ch_ref_genome,
            )
            VcfToMatrixTable(
                channel.fromPath(params.merged_vcf, checkIfExists: true),
                ch_merged_bed,
            )
            ch_mt = VcfToMatrixTable.out
        }
	}

	// new stuff from here - we have a MT, now export sites only fragments from it
	MakeSitesOnlyVcfsFromMt(ch_mt)

	// MakeSitesOnlyVcfsFromMt.out will be a channel of *.bgz vcf fragments

    // apply the gnomAD and AlphaMissense annotations using echtvar
    AnnotateWithEchtvar(
        MakeSitesOnlyVcfsFromMt.out.collect().flatten(),
        ch_gnomad_zip,
        ch_alphamissense_zip,
    )

    // annotate transcript consequences with bcftools csq
    AnnotateCsqWithBcftools(
        AnnotateWithEchtvar.out.collect().flatten(),
        ch_gff,
        ch_ref_genome,
    )

    // reformat the annotations in the VCF, retain as a Hail Table
    ReformatAnnotatedVcfIntoHailTable(
        AnnotateCsqWithBcftools.out.collect().flatten(),
        ch_bed,
        ch_mane,
    )

    // combine the join-VCF and annotations as a HailTable
    TransferAnnotationsToMatrixTable(
        ch_mt,
        ReformatAnnotatedVcfIntoHailTable.out.collect(),
    )
}
