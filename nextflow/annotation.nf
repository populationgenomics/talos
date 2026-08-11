#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
This workflow is the annotation process for the Talos pipeline.

It requires a single VCF, which can be single- or multi-sample. This is then annotated and reformatted for Talos.
The specific annotations are:

- gnomAD v4.1 frequencies and alphamissense annotations, applied to the joint VCF using echtvar
- Transcript consequences, using BCFtools annotate

Everything here is cohort-agnostic and cacheable. Output is the annotated VCF shards themselves -
MANE and ClinVar are applied in the per-run workflow, where the labelling happens.
*/

include { AnnotateCsqWithBcftools } from './modules/annotation/AnnotateCsqWithBcftools/main'
include { AnnotateWithEchtvar } from './modules/annotation/AnnotateWithEchtvar/main'
include { MergeVcfsWithBcftools } from './modules/annotation/MergeVcfsWithBcftools/main'
include { NormaliseVcf } from './modules/annotation/NormaliseVcf/main'
include { SplitVcf } from './modules/annotation/SplitVcf/main'


workflow ANNOTATION {
	take:
		ch_gff
		ch_ref_genome
		ch_inputs

    main:
    // populate various input channels - these are downloaded by the large_files/gather_files.sh script, or the prep wf
    if (!file(params.alphamissense_zip).exists()) {
        println "AlphaMissense data must be encoded for echtvar, run the Talos Prep workflow (talos_preparation.nf)"
        exit 1
    }
    ch_alphamissense_zip = channel.fromPath(params.alphamissense_zip, checkIfExists: true).first()

    ch_gnomad_zip = channel.fromPath(params.gnomad_zip, checkIfExists: true).first()

    ch_inputs_branched = ch_inputs.branch {
        shards: it[2] == 'shards'
        vcf_dir: it[2] == 'ss_vcf_dir'
        single_vcf: it[2] == 'vcf'
    }

    // Process shards
    ch_from_shards = ch_inputs_branched.shards.flatMap { cohort, path, _type ->
        def vcfs = files("${path}/*.${params.input_vcf_extension}")
        vcfs.collect { vcf -> tuple(cohort, vcf) }
    }

    // Process single-sample components
    ch_vcf_dir_inputs = ch_inputs_branched.vcf_dir.map { cohort, path, _type ->
        def vcfs = files("${path}/*.${params.input_vcf_extension}")
        def tbis = vcfs.collect { file("${it}.tbi") }
        tuple(cohort, vcfs, tbis)
    }

    MergeVcfsWithBcftools(
        ch_vcf_dir_inputs,
        ch_ref_genome,
    )
    ch_merged_vcfs = MergeVcfsWithBcftools.out

    // Process single VCF
    ch_single_vcfs = ch_inputs_branched.single_vcf.map { cohort, path, _type ->
        tuple(cohort, file(path, checkIfExists: true))
    }

    // Combine single VCFs and merged VCFs
    ch_to_split = ch_single_vcfs.mix(ch_merged_vcfs)

    if ((params.vcf_split_n ?: 0) > 0) {
        SplitVcf(ch_to_split)
        ch_vcfs = SplitVcf.out.flatMap { cohort, splits ->
            def split_list = splits instanceof Collection ? splits : [splits]
            split_list.collect { tuple(cohort, it) }
        }
    } else {
        ch_vcfs = ch_to_split
    }

    // mix in shards
    ch_all_vcfs = ch_vcfs.mix(ch_from_shards)

	NormaliseVcf(
        ch_all_vcfs,
        ch_ref_genome,
    )

	AnnotateWithEchtvar(
        NormaliseVcf.out,
        ch_gnomad_zip,
        ch_alphamissense_zip,
    )

    // annotate transcript consequences with bcftools csq
    AnnotateCsqWithBcftools(
        AnnotateWithEchtvar.out,
        ch_gff,
        ch_ref_genome,
    )

    emit:
        // one entry per shard, as [cohort, vcf] - the per-run workflow keeps the scatter open
    	vcfs = AnnotateCsqWithBcftools.out
}
