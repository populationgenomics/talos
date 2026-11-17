process AnnotateShard {
    container params.container
    publishDir { annotated_dir }, mode: 'copy'

    // annotate one shard. echtvar (gnomAD + AlphaMissense) and bcftools csq run on the sites-only
    // VCF, which is a small fraction of the full-width shard for any multi-sample cohort. The INFO
    // fields they add are then lifted onto the full-width BCF with bcftools annotate, matching on
    // CHROM/POS/REF/ALT - exact, as NormaliseVcf has already atomised every record.
    //
    // The sites file is unfiltered, so every full-width record has a match and the lift is a plain
    // INFO copy. The common-variant cut (gnomAD AF >= 0.05, the loosest threshold Talos applies)
    // is made on the full-width stream afterwards - echtvar writes gnomad_AF=0 for variants absent
    // from gnomAD, and -e also retains a record if the field were ever missing, so novel variants
    // are always kept.
    //
    // csq's per-sample FORMAT/BCSQ bitmask is not produced (there are no samples in the sites file)
    // - nothing downstream reads it, only INFO/BCSQ
    input:
        tuple val(cohort), path(full_bcf), path(full_csi), path(sites_vcf), val(annotated_dir)
        path gnomad_zip
        path am_zip
        path gff3
        path reference

    output:
        tuple val(cohort), path("${full_bcf.simpleName}_csq.vcf.bgz"), path("${full_bcf.simpleName}_csq.vcf.bgz.tbi")

    script:
        def name = full_bcf.simpleName
        """
        set -euo pipefail

        echtvar anno \
            -e ${gnomad_zip} \
            -e ${am_zip} \
            ${sites_vcf} \
            "${name}_sites_echtvar.vcf.bgz"

        bcftools csq --force -f "${reference}" \
            --greedy 1 \
            --local-csq \
            -g ${gff3} \
            --unify-chr-names 'chr,-,chr' \
            -B 20 \
            -Oz -o "${name}_sites_csq.vcf.bgz" \
            "${name}_sites_echtvar.vcf.bgz"
        tabix "${name}_sites_csq.vcf.bgz"

        # -c INFO lifts every INFO field (and its header line) from the annotated sites onto the
        # full-width records, then drop common variants. This output is published and re-read every
        # reanalysis cycle, so it keeps the default compression level - but compresses in parallel
        bcftools annotate \
            -a "${name}_sites_csq.vcf.bgz" \
            -c INFO \
            -Ou \
            --no-version \
            ${full_bcf} | \
        bcftools view \
            -e 'INFO/gnomad_AF >= 0.05' \
            -Oz \
            --threads ${task.cpus} \
            --no-version \
            -o "${name}_csq.vcf.bgz" \
            -W=tbi
        """
}
