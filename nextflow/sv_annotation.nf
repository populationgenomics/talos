#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
This workflow is the Structural Variant annotation process for the Talos pipeline.

It takes a joint-called SV VCF per cohort, supplied as the `sv` column of the input TSV, and annotates it for
Talos. The specific annotations are:

- gene consequences, using GATK SVAnnotate against the MANE GTF. PREDICTED_LOF is the field Talos requires
- gnomAD v4.1 population frequencies, using SVAFotate against its published popAF BED
- a rename of SVAFotate's field names to the ones talos.run_sv_filtering reads

This is deliberately separate from the ANNOTATION workflow rather than a branch inside it. That workflow's
normalise, region-filter, split and echtvar steps are all SNV-specific and would each need bypassing, and the
SV input is a single joint-called VCF with no equivalent of the shards / ss_vcf_dir branching.

Every step is gated per-cohort. A cohort is skipped entirely if it has no SV VCF, or if its final annotated
VCF is already present in the cohort's annotated directory - SVAnnotate and SVAFotate are both expensive enough, and the
inputs stable enough, that re-running them on an unchanged callset is pure waste.
*/

include { AnnotateSvWithGatk } from './modules/annotation/AnnotateSvWithGatk/main'
include { AnnotateSvWithSvafotate } from './modules/annotation/AnnotateSvWithSvafotate/main'
include { CreateSequenceDictionary } from './modules/prep/CreateSequenceDictionary/main'
include { RenameSvAfFields } from './modules/annotation/RenameSvAfFields/main'
include { PrepareSvForSvAnnotate } from './modules/annotation/PrepareSvForSvAnnotate/main'


workflow SV_ANNOTATION {
    take:
        ch_ref_genome
        // one entry per cohort: [cohort, sv_vcf, config, pedigree, annotated_dir] - the VCF is cut down to the
        // pedigree's samples before annotation
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

    // previously annotated SV VCFs are read from the cohort's annotated directory, alongside its small-variant shards.
    // keep this a Path and resolve() against it - interpolating a cloud Path into a string
    // drops the gs:// scheme, so the exists() check silently ran against the local filesystem

    // split cohorts three ways - no SV data ([] is falsy), already annotated, or work to do
    // row is tuple(cohort, sv_vcf, config, pedigree, annotated_dir)
    ch_branched = ch_sv_inputs.branch { row ->
        no_sv:    !row[1]
        complete: file(row[4]).resolve("${row[0]}_sv_annotated.vcf.bgz").exists()
        pending:  true
    }

    // pick the previously annotated VCF up from the annotated directory instead of regenerating it
    ch_complete = ch_branched.complete.map { row ->
        def annotated = file(row[4]).resolve("${row[0]}_sv_annotated.vcf.bgz")
        println "Annotated SV VCF for ${row[0]} already exists (${annotated.toUriString()}), skipping SV annotation"
        def tbi = annotated.resolveSibling("${annotated.name}.tbi")
        if (!tbi.exists()) {
            error("Annotated SV VCF for ${row[0]} exists but its index is missing: ${tbi.toUriString()}")
        }
        tuple(row[0], annotated, tbi)
    }

    ch_pending_vcfs = ch_branched.pending.map { row ->
        tuple(row[0], row[1], file("${row[1]}.tbi", checkIfExists: true), row[3])
    }

    // the config and annotated dir travel with each cohort, and are only needed for the (published) rename step
    ch_pending_configs = ch_branched.pending.map { row -> tuple(row[0], row[2], row[4]) }

    // SVAnnotate needs a sequence dictionary, which large_files does not ship alongside ref.fa.
    // Building it is gated on there being pending work - `first()` over an empty channel emits nothing, so
    // when every cohort was skipped CreateSequenceDictionary receives no input and never runs
    ch_work_pending = ch_pending_vcfs.first()

    if (!file(params.ref_dict).exists()) {
        CreateSequenceDictionary(
            ch_ref_genome.combine(ch_work_pending).map { ref_fa, _cohort, _vcf, _tbi, _pedigree -> ref_fa },
        )
        ch_ref_dict = CreateSequenceDictionary.out.first()
    } else {
        ch_ref_dict = channel.fromPath(params.ref_dict, checkIfExists: true).first()
    }
    // SVAnnotate aborts the entire run on the first complex variant whose CPX_INTERVALS are not in
    // coordinate order, which GATK-SV's delINVdel records never are - so sort them before GATK sees them.
    // SVAnnotate also only works with symbolic alleles, so the same pass removes any non-symbolic ALT record
    PrepareSvForSvAnnotate(ch_pending_vcfs)
    AnnotateSvWithGatk(
        PrepareSvForSvAnnotate.out,
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
        // every cohort with SV data, newly annotated and reused alike
        annotated = RenameSvAfFields.out.mix(ch_complete)
}
