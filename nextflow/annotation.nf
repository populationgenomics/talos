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

include { SplitVcf } from './modules/annotation/SplitVcf/main'
include { AnnotateCsqWithBcftools } from './modules/annotation/AnnotateCsqWithBcftools/main'
include { AnnotateWithEchtvar } from './modules/annotation/AnnotateWithEchtvar/main'
include { AnnotatedVcfIntoMatrixTable } from './modules/annotation/AnnotatedVcfIntoMatrixTable/main'
include { CreateRoiFromGff3 } from './modules/annotation/CreateRoiFromGff3/main'
include { EncodeAlphaMissense } from './modules/annotation/EncodeAlphaMissense/main'
include { NormaliseAndRegionFilterVcf } from './modules/annotation/NormaliseAndRegionFilterVcf/main'
include { ParseAlphaMissense } from './modules/annotation/ParseAlphaMissense/main'
include { ParseManeIntoJson } from './modules/annotation/ParseManeIntoJson/main'

workflow {

    // pull and parse the MANE data into JSON file
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

    // populate ref genome input channel
    ch_ref_genome = channel.fromPath(
    	params.ref_genome,
    	checkIfExists: true,
	)

    // generate the AlphaMissense Echtvar annotation source - long running, but only requires generation once
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

    // generate the Region-of-interest BED file from Ensembl GFF3
    // generates a per-gene BED file with ID annotations
    // and a overlap-merged version of the same for more efficient region filtering
    // use of this BED file limits the region of interest to genic areas. To be reconsidered if we expand to non-genic.
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

	// require a pre-merged callset
	ch_merged_vcf = channel.fromPath(params.merged_vcf, checkIfExists: true)
	SplitVcf(ch_merged_vcf)

	NormaliseAndRegionFilterVcf(
        SplitVcf.out.flatten(),
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

    // reformat the annotations in the VCF, retain as a Hail Table
    AnnotatedVcfIntoMatrixTable(
        AnnotateCsqWithBcftools.out,
        ch_bed.first(),
        ch_mane.first(),
    )
}
