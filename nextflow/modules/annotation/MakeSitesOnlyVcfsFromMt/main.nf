
process MakeSitesOnlyVcfsFromMt {
    container params.container

    // take the merged VCF and index
    input:
        path mt

    // strip this VCF down to sites-only
    publishDir params.cohort_output_dir

    output:
        path "*.bgz"

    script:
    """
    python -m talos.annotation_scripts.MakeSitesOnlyVcfsFromMt \
        --input ${mt} \
        --output here.vcf.bgz
    mv here.vcf.bgz/part*.bgz .
    rm -r here.vcf.bgz
    """
}
