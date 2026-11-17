process MakeScatterRegions {
    container params.container

    // derive scatter regions from the VCF index alone - bcftools index -s reads contig, length
    // and record count straight from the tbi, so the VCF itself is never staged. Each contig is
    // cut into ceil(records / vcf_split_n) equal-length windows, keeping shard record counts near
    // the target. Contigs with no length available (bcftools reports ".") are emitted whole.
    // Only canonical contigs (chr1-22, X, Y) are scattered; alts, decoys, and unplaced/unlocalised
    // contigs are dropped.
    // A separate Mito pathway exists ('mito' in the input file) with Mito-genome sensitive consequence annotation
    // so chrM is also removed at this point.
    input:
        tuple val(cohort), path(tbi)

    output:
        tuple val(cohort), path("${cohort}_scatter_regions.txt")

    script:
        """
        set -euo pipefail

        bcftools index -s ${tbi} | \
        awk -v target=${params.vcf_split_n} '
            {
                chr = \$1; len = \$2; n = \$3
                if (n == 0) next
                if (chr !~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)\$/) next
                k = int((n + target - 1) / target)
                if (k <= 1 || len == ".") { print chr; next }
                step = int((len + k - 1) / k)
                for (start = 1; start <= len; start += step) {
                    end = start + step - 1
                    if (end > len) end = len
                    printf "%s:%d-%d\\n", chr, start, end
                }
            }' > "${cohort}_scatter_regions.txt"

        # an empty regions file means an empty (or unindexable) VCF - fail here, not downstream
        [ -s "${cohort}_scatter_regions.txt" ]
        """
}
