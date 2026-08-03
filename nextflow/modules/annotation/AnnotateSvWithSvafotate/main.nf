process AnnotateSvWithSvafotate {
    container params.svafotate_container

    input:
        tuple val(cohort), path(vcf), path(tbi)
        path svafotate_bed

    // writes an uncompressed VCF - this image carries no htslib CLI, so compression and indexing
    // are left to RenameSvAfFields, which runs in the Talos container
    output:
        tuple val(cohort), path("${cohort}_sv_popaf.vcf")

    script:
        """
        set -euo pipefail

        # small interjection here - the svafotate process requires a HUGE BED file, with multiple sources and
        # ~170 columns. This is memory intensive enough that I cannot complete it locally... this pre-processing
        # reduces the resources required for annotation
        #
        # single streaming pass: decompress, keep the header plus gnomAD rows, recompress. Do NOT use
        # `zcat ... | head -1` for the header - head exits after one line, zcat takes SIGPIPE on the rest of a
        # multi-GB file, and `set -o pipefail` turns that into exit 141. `NR == 1 ||` also stops the header
        # being emitted twice: its column names contain "gnomAD", so a bare grep matches it as a data row.
        #
        # \\\$6 is the SOURCE column. Nextflow script blocks are Groovy GStrings, so a bare \$6 is interpolated away
        # to an empty string and awk dies with a syntax error on `NR == 1 ||  == "gnomAD"`.
        # -F is set for the same reason: awk's default separator collapses runs of whitespace, so a future BED with an
        # empty field before SOURCE would shift the numbering silently.
        zcat ${svafotate_bed} \
            | awk -F'\\t' 'NR == 1 || \$6 == "gnomAD"' \
            | gzip -c > new_af_file.bed.gz

        # -s gnomAD is required. Without it Max_AF is the maximum across the BED's four merged callsets
        # (gnomAD, CCDG, TOPMed, ThousG), which would be mislabelled as a gnomAD frequency downstream.
        #
        # -f must be set explicitly. It defaults to 0.001, which is effectively "any single-basepair overlap
        # matches", and would annotate a 200bp deletion with a 3Mb gnomAD deletion's frequency.
        #
        # The reference BED is passed exactly as downloaded. SVAFotate strips the `chr` prefix from the query
        # VCF but not from this BED, so its Ensembl-style contig names are already correct - adding a `chr`
        # prefix here would silently zero every frequency in the callset.
        svafotate annotate \
            -v ${vcf} \
            -o ${cohort}_sv_popaf.vcf \
            -b new_af_file.bed.gz \
            -s gnomAD \
            -f ${params.sv_overlap_fraction} \
            -a best \
            --ins \
            --cpu ${task.cpus}
        """
}
