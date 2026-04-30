#!/usr/bin/env nextflow

/*
This workflow takes ownership over all the input data preparation steps across both Talos workflows (annotation + talos)

The data prepared in this workflow is:

- AlphaMissense data, reformatted into a HailTable
- PanelApp data, freshly downloaded
- ClinVar data, freshly downloaded and reformatted into Hail Tables
- Mane data, downloaded and reformatted
- BED file representing ROI for all Ensembl genes
*/

nextflow.enable.dsl=2

include { AnnotateClinvarWithBcftools } from './nextflow/modules/prep/AnnotateClinvarWithBcftools/main'
include { CreateRoiFromGff3 } from './nextflow/modules/prep/CreateRoiFromGff3/main'
include { DownloadClinVarFiles } from './nextflow/modules/prep/DownloadClinVarFiles/main'
include { DownloadPanelApp } from './nextflow/modules/prep/DownloadPanelApp/main'
include { EncodeAlphaMissense } from './nextflow/modules/prep/EncodeAlphaMissense/main'
include { EncodeMitimpact } from './nextflow/modules/prep/EncodeMitimpact/main'
include { EncodeMitotip } from './nextflow/modules/prep/EncodeMitotip/main'
include { EncodeNapogee } from './nextflow/modules/prep/EncodeNapogee/main'
include { MakeClinvarbitrationPm5 } from './nextflow/modules/prep/MakeClinvarbitrationPm5/main'
include { ResummariseRawSubmissions } from './nextflow/modules/prep/ResummariseRawSubmissions/main'
include { ParseManeIntoJson } from './nextflow/modules/prep/ParseManeIntoJson/main'


workflow {
    main:
    def timestamp = new java.util.Date().format('yyyy-MM')

    ch_gff = channel.fromPath(params.ensembl_gff, checkIfExists: true)

    // generate the AlphaMissense HT - long running, stored in a separate folder
    if (!file(params.alphamissense_zip).exists()) {
    	ch_alphamissense_tsv = channel.fromPath(params.alphamissense_tsv, checkIfExists: true)
        EncodeAlphaMissense(ch_alphamissense_tsv)
        ch_alphamissense_zip = EncodeAlphaMissense.out
    } else {
        ch_alphamissense_zip = channel.fromPath(params.alphamissense_zip)
    }

    // generate/format the various Mitochondrial annotation sources
    if (!file(params.mitimpact_zip).exists()) {
    	ch_mitimpact_tsv = channel.fromPath(params.mitimpact_tsv, checkIfExists: true)
        EncodeMitimpact(ch_mitimpact_tsv)
        ch_mitimpact_zip = EncodeMitimpact.out
    } else {
        ch_mitimpact_zip = channel.fromPath(params.mitimpact_zip)
    }

    // generate/format the various Mitochondrial annotation sources
    if (!file(params.mitotip_zip).exists()) {
    	ch_mitotip_tsv = channel.fromPath(params.mitotip_tsv, checkIfExists: true)
        EncodeMitotip(ch_mitotip_tsv)
        ch_mitotip_zip = EncodeMitotip.out
    } else {
        ch_mitotip_zip = channel.fromPath(params.mitotip_zip)
    }

    // generate/format the various Mitochondrial annotation sources
    if (!file(params.napogee_zip).exists()) {
    	ch_napogee_tsv = channel.fromPath(params.napogee_tsv, checkIfExists: true)
        EncodeNapogee(ch_napogee_tsv)
        ch_napogee_zip = EncodeNapogee.out
    } else {
        ch_napogee_zip = channel.fromPath(params.napogee_zip)
    }

    // does this month's clinvarbitration data exist?
    String current_clinvarbitration_all = "${params.processed_annotations}/clinvarbitration_${timestamp}.ht"
    String current_clinvarbitration_pm5 = "${params.processed_annotations}/clinvarbitration_${timestamp}.pm5.ht"

    if (file(current_clinvarbitration_pm5).exists()) {
        ch_clinvar_all = channel.fromPath(current_clinvarbitration_all)
        ch_clinvar_pm5 = channel.fromPath(current_clinvarbitration_pm5)
    } else {
        // new workflow elements to go and create it from raw data
        String subfile = "${params.large_files}/submissions_${timestamp}.txt.gz"
        String varfile = "${params.large_files}/variants_${timestamp}.txt.gz"

        if (file(subfile).exists() && file(varfile).exists()) {
            ch_clinvar_sub = channel.fromPath(subfile)
            ch_clinvar_var = channel.fromPath(varfile)
        } else {
            println "Attempting to download ClinVar raw data, requires internet connection."
            println "If this step fails, try re-running gather_files.sh in the `large_files` directory."
            // this step requires an internet connection, which may be problematic at some sites
            DownloadClinVarFiles(timestamp)
            ch_clinvar_sub = DownloadClinVarFiles.out.submissions
            ch_clinvar_var = DownloadClinVarFiles.out.variants
        }

        // reinterpret the results using altered heuristics
        ResummariseRawSubmissions(
            ch_clinvar_var,
            ch_clinvar_sub,
            timestamp,
        )

        ch_ref_fa = channel.fromPath(params.ref_genome, checkIfExists: true)

        // annotate the SNV VCF using BCFtools
        AnnotateClinvarWithBcftools(
            ResummariseRawSubmissions.out.vcf,
            ch_ref_fa,
            ch_gff,
        )

        MakeClinvarbitrationPm5(
            AnnotateClinvarWithBcftools.out,
            timestamp,
        )

        ch_clinvar_all = ResummariseRawSubmissions.out.ht
        ch_clinvar_pm5 = MakeClinvarbitrationPm5.out
    }

    // generate the Region-of-interest BED file from Ensembl GFF3
    // generates a per-gene BED file with ID annotations
    // and a overlap-merged version of the same for more efficient region filtering
    if (!file(params.ensembl_bed).exists() || !file(params.ensembl_merged_bed).exists()) {
        CreateRoiFromGff3(ch_gff)
        ch_bed = CreateRoiFromGff3.out.bed
        ch_merged_bed = CreateRoiFromGff3.out.merged_bed
    } else {
        ch_merged_bed = channel.fromPath(params.ensembl_merged_bed, checkIfExists: true)
        ch_bed = channel.fromPath(params.ensembl_bed, checkIfExists: true)
    }

    // pull and parse the MANE data into a Hail Table
    if (!file(params.mane_json).exists()) {
        ch_mane_summary = channel.fromPath(params.mane, checkIfExists: true)
        ParseManeIntoJson(ch_mane_summary)
        ch_mane_json = ParseManeIntoJson.out
    } else {
        ch_mane_json = channel.fromPath(params.mane_json, checkIfExists: true)
    }

    String current_panelapp = "${params.processed_annotations}/panelapp_${timestamp}.json"

    if (!file(current_panelapp).exists()) {
        DownloadPanelApp(
            ch_mane_json,
            timestamp,
        )
        panelapp_out = DownloadPanelApp.out
    } else {
        panelapp_out = channel.fromPath(current_panelapp, checkIfExists: true)
    }

    // use workflow outputs, not individual copies
    publish:
        alphamissense = ch_alphamissense_zip
        bed = ch_bed
        merged_bed = ch_merged_bed
        clinvar_all = ch_clinvar_all
        clinvar_pm5 = ch_clinvar_pm5
        mitimpact = ch_mitimpact_zip
        mitotip = ch_mitotip_zip
        mane_json = ch_mane_json
        napogee = ch_napogee_zip
        panelapp_out = panelapp_out
}

output {
    alphamissense {
    }
    bed {
    }
    merged_bed {
    }
    clinvar_all {
    }
    clinvar_pm5 {
    }
    mane_json {
    }
    mitimpact {
    }
    mitotip {
    }
    napogee {
    }
    panelapp_out {
    }
}
