#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
This workflow is the Structural Variant annotation process for the Talos pipeline.

It takes a joint-called SV VCF per cohort, supplied as the `sv` column of the input TSV, and annotates it for
Talos. The specific annotations are:

- gene consequences, using GATK SVAnnotate against the MANE GTF. PREDICTED_LOF is the field Talos requires
- gnomAD v4.1 population frequencies, using SVAFotate against its published popAF BED
- a rename of SVAFotate's field names to the ones talos.run_hail_filtering_sv reads

This is deliberately separate from the ANNOTATION workflow rather than a branch inside it. That workflow's
normalise, region-filter, split and echtvar steps are all SNV-specific and would each need bypassing, and the
SV input is a single joint-called VCF with no equivalent of the shards / ss_vcf_dir branching.

Every step is gated per-cohort. A cohort is skipped entirely if it has no SV VCF, or if its final annotated
VCF is already present in the output directory - SVAnnotate and SVAFotate are both expensive enough, and the
inputs stable enough, that re-running them on an unchanged callset is pure waste.
*/

include { AnnotateSvWithGatk } from './modules/annotation/AnnotateSvWithGatk/main'
include { AnnotateSvWithSvafotate } from './modules/annotation/AnnotateSvWithSvafotate/main'
include { CreateSequenceDictionary } from './modules/prep/CreateSequenceDictionary/main'
include { RenameSvAfFields } from './modules/annotation/RenameSvAfFields/main'
include { SortCpxIntervals } from './modules/annotation/SortCpxIntervals/main'


workflow SV_ANNOTATION {
    take:
        ch_ref_genome
        ch_sv_inputs

    main:
    // both of these are downloaded by large_files/gather_files.sh
    if (!file(params.mane_gtf).exists()) {
        println "MANE GTF (${params.mane_gtf}) is not available, please re-run the input file download script"
        exit 1
    }
    ch_mane_gtf = channel.fromPath(params.mane_gtf, checkIfExists: true).first()

    if (!file(params.svannotate_noncoding_bed).exists()) {
        println "Non-coding BED (${params.svannotate_noncoding_bed}) is not available, please re-run the input file download script"
        exit 1
    }
    ch_noncoding_bed = channel.fromPath(params.svannotate_noncoding_bed, checkIfExists: true).first()

    if (!file(params.svafotate_bed).exists()) {
        println "SVAFotate BED (${params.svafotate_bed}) is not available, please re-run the input file download script"
        exit 1
    }
    ch_svafotate_bed = channel.fromPath(params.svafotate_bed, checkIfExists: true).first()

    // split cohorts three ways - no SV data, already annotated, or work to do
    // row is tuple(cohort, sv_vcf, config)
    ch_branched = ch_sv_inputs.branch { row ->
        sentinel: row[1].name == 'NO_SV'
        complete: file("${workflow.outputDir}/${row[0]}_outputs/${row[0]}_sv_annotated.vcf.bgz").exists()
        pending:  true
    }

    // pick the previously annotated VCF up from the output directory instead of regenerating it
    ch_complete = ch_branched.complete.map { row ->
        def annotated = "${workflow.outputDir}/${row[0]}_outputs/${row[0]}_sv_annotated.vcf.bgz"
        println "Annotated SV VCF for ${row[0]} already exists (${annotated}), skipping SV annotation"
        tuple(row[0], file(annotated), file("${annotated}.tbi", checkIfExists: true))
    }

    ch_pending_vcfs = ch_branched.pending.map { row ->
        tuple(row[0], row[1], file("${row[1]}.tbi", checkIfExists: true))
    }

    // the config travels with each cohort, and is only needed for the rename step
    ch_pending_configs = ch_branched.pending.map { row -> tuple(row[0], row[2]) }

    // SVAnnotate needs a sequence dictionary, which large_files does not ship alongside ref.fa.
    // Building it is gated on there being pending work - `first()` over an empty channel emits nothing, so
    // when every cohort was skipped CreateSequenceDictionary receives no input and never runs
    ch_work_pending = ch_pending_vcfs.first()

    if (!file(params.ref_dict).exists()) {
        CreateSequenceDictionary(
            ch_ref_genome.combine(ch_work_pending).map { ref_fa, _cohort, _vcf, _tbi -> ref_fa },
        )
        ch_ref_dict = CreateSequenceDictionary.out.first()
    } else {
        ch_ref_dict = channel.fromPath(params.ref_dict, checkIfExists: true).first()
    }
    // SVAnnotate aborts the entire run on the first complex variant whose CPX_INTERVALS are not in
    // coordinate order, which GATK-SV's delINVdel records never are - so sort them before GATK sees them
    SortCpxIntervals(ch_pending_vcfs)
    AnnotateSvWithGatk(
        SortCpxIntervals.out,
        ch_mane_gtf,
        ch_noncoding_bed,
        ch_ref_dict,
    )

    AnnotateSvWithSvafotate(
        AnnotateSvWithGatk.out,
        ch_svafotate_bed,
    )

    RenameSvAfFields(
        AnnotateSvWithSvafotate.out.join(ch_pending_configs),
    )

    emit:
        annotated = RenameSvAfFields.out.mix(ch_complete)
}
