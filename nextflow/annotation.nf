#!/usr/bin/env nextflow

/*
This workflow is the annotation process for the Talos pipeline.

It requires a single VCF, which can be single- or multi-sample. This is then annotated and reformatted for Talos.
The specific annotations are:

- gnomAD v4.1 frequencies and alphamissense annotations, applied to the joint VCF using echtvar
- Transcript consequences, using BCFtools annotate
- MANE trancript IDs and corresponding ENSP IDs, applied using Hail
*/

nextflow.enable.dsl=2

include { AnnotateCsqWithBcftools } from './modules/annotation/AnnotateCsqWithBcftools/main'
include { AnnotatedVcfIntoMatrixTable } from './modules/annotation/AnnotatedVcfIntoMatrixTable/main'
include { AnnotateWithEchtvar } from './modules/annotation/AnnotateWithEchtvar/main'
include { MergeVcfsWithBcftools } from './modules/annotation/MergeVcfsWithBcftools/main'
include { NormaliseAndRegionFilterVcf } from './modules/annotation/NormaliseAndRegionFilterVcf/main'
include { SplitVcf } from './modules/annotation/SplitVcf/main'


workflow ANNOTATION {
	take:
		ch_gff
		ch_mane
		ch_ref_genome
		ch_inputs

    main:
    // populate various input channels - these are downloaded by the large_files/gather_file.sh script, or the prep wf
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

    ch_gnomad_zip = channel.fromPath(params.gnomad_zip, checkIfExists: true)

    ch_inputs_branched = ch_inputs.branch {
        shards: it[2] == 'shards'
        vcf_dir: it[2] == 'ss_vcf_dir'
        single_vcf: it[2] == 'vcf'
    }

    // Process shards
    ch_from_shards = ch_inputs_branched.shards.flatMap { cohort, path, type ->
        def vcfs = file("${path}/*.${params.input_vcf_extension}")
        vcfs.collect { vcf -> tuple(cohort, vcf) }
    }

    // Process single-sample components
    ch_vcf_dir_inputs = ch_inputs_branched.vcf_dir.map { cohort, path, type ->
        def vcfs = file("${path}/*.${params.input_vcf_extension}")
        def tbis = vcfs.collect { file("${it}.tbi") }
        tuple(cohort, vcfs, tbis)
    }

    MergeVcfsWithBcftools(
        ch_vcf_dir_inputs,
        ch_ref_genome.first(),
    )
    ch_merged_vcfs = MergeVcfsWithBcftools.out.merged

    // Process single VCF
    ch_single_vcfs = ch_inputs_branched.single_vcf.map { cohort, path, type ->
        tuple(cohort, file(path, checkIfExists: true))
    }

    // Combine single VCFs and merged VCFs
    ch_to_split = ch_single_vcfs.mix(ch_merged_vcfs)

    // decide whether to split and parallelise, or run as a single operation
    // if config value is absent completely, skip this step
    if ((params.vcf_split_n ?: 0) > 0) {
        SplitVcf(ch_to_split)
        ch_vcfs = SplitVcf.out.flatMap { cohort, splits -> splits.collect { tuple(cohort, it) } }
    } else {
        ch_vcfs = ch_to_split
    }

    // mix in shards
    ch_all_vcfs = ch_vcfs.mix(ch_from_shards)

	NormaliseAndRegionFilterVcf(
        ch_all_vcfs,
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

    emit:
    	mts = AnnotatedVcfIntoMatrixTable.out.groupTuple(by: 0)
}
