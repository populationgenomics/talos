process ConcatLabelledVcfs {
    container params.container
    publishDir { analysis_dir }, mode: 'copy'

    // gather the labelled shards back into the single per-cohort VCF ValidateMOI expects.
    // each shard is coordinate-sorted internally, but shard names don't sort genomically
    // (region tags sort lexically - chr10 before chr2 - and pre-sharded inputs are named by
    // whoever made them), so the concatenated stream is re-sorted rather than trusted
    input:
        // indexes aren't named in the script, but staging them is what lets `concat -a` work
        tuple val(cohort), path(vcfs), path(indexes), val(analysis_dir)

    output:
        tuple val(cohort), path("${cohort}_small_variants_labelled.vcf.bgz"), path("${cohort}_small_variants_labelled.vcf.bgz.tbi")

    script:
        def vcf_list = (vcfs instanceof Collection ? vcfs.sort{ vcf -> vcf.name } : [vcfs]).join(' ')
        """
        set -euo pipefail

        mkdir -p sort_tmp

        bcftools concat \
            -a \
            --no-version \
            -Ou \
            ${vcf_list} | \
        bcftools sort \
            -T sort_tmp \
            -Oz \
            -o "${cohort}_small_variants_labelled.vcf.bgz" \
            -W=tbi -
        """
}
