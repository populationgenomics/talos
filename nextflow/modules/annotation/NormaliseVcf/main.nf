process NormaliseVcf {
    container params.container

    // normalise one scatter region of the cohort VCF - this is where the scatter actually happens,
    // reading only the region's blocks through the index. --regions-overlap 0 assigns a record to
    // the shard containing its POS, so a deletion spanning a region boundary is emitted exactly once.
    // region 'all' streams the whole file (pre-sharded inputs, or scatter disabled) and needs no
    // index - tbi is [] on that path, staging nothing.
    // ref_genome here is used to create parsimonious representations
    //
    // Two outputs per shard: the full-width normalised BCF (all samples, only ever read again by the
    // lift-over in AnnotateShard - BCF because it is cheaper to write and re-read than bgzipped text;
    // indexed because bcftools annotate needs an index on both sides when -a is a VCF),
    // and a sites-only VCF (no FORMAT columns) which is what echtvar and bcftools csq actually run on.
    // Both come from the same normalised stream, so records correspond 1:1 by CHROM/POS/REF/ALT.
    // fill-tags runs before the sites are extracted, as AC/AF/AN need the genotypes. It is skipped when
    // params.recalculate_ac_af_an is false, keeping AC/AF/AN already in the VCF (e.g. from a larger callset).
    //
    // The region is cut down to the cohort's pedigree samples before anything else reads it, dropping any
    // allele (-a) and record (-c 1) no remaining sample carries - a dataset-wide callset is annotated, stored
    // and filtered at cohort size. A VCF already holding only pedigree samples is read as-is. Pedigree samples
    // absent from the VCF are a warning (the Python stages tolerate them too), no shared samples at all an error.
    // When recalculating, AC/AF/AN are written by fill-tags below, so describe the cohort, not the source
    // callset. When not, the subset passes -I so bcftools view doesn't overwrite INFO/AC,AN with the
    // subset's counts - -a and -c 1 still work from the subset genotypes, and trim Number=A fields to match.
    //
    // FORMAT/PL is dropped here, at the head of the pipeline: Talos never reads it (GT, GQ, DP, AD
    // and the phasing fields are all it uses) and it is the widest per-sample field, so removing it
    // shrinks every downstream pass and the published shards by roughly a third.
    // The full-width BCF is written once and read once, so it gets the lightest deflate level (-Ob1)
    // with threaded compression - deflate is most of this task's runtime, not norm itself
    input:
        tuple val(cohort), path(vcf), path(tbi), val(region), path(pedigree)
        path ref_genome

    output:
        tuple val(cohort), path("*_normalised.bcf"), path("*_normalised.bcf.csi"), path("*_sites.vcf.bgz")

    script:
        def out_name = region == 'all' ? vcf.simpleName : "${vcf.simpleName}_${region.replaceAll(/[:\-]/, '_')}"
        // --regions-overlap 0 = POS-in-region only, so boundary-spanning records land in exactly one shard
        def region_args = region == 'all' ? '' : "--regions ${region} --regions-overlap 0"
        def norm_args = "-m -any -f ${ref_genome} -Ou --no-version"
        def recalculate = params.recalculate_ac_af_an
        // without recalculation, -I stops the sample subset rewriting INFO/AC and INFO/AN
        def subset_args = recalculate ? '' : '-I'
        // fill-tags or a plain conversion - both write the indexed full-width BCF
        def write_bcf = recalculate
            ? "bcftools +fill-tags -Ob1 --threads ${task.cpus} --no-version -o \"${out_name}_normalised.bcf\" -W - -- -t AC,AF,AN"
            : "bcftools view -Ob1 --threads ${task.cpus} --no-version -o \"${out_name}_normalised.bcf\" -W -"
        """
        set -euo pipefail

        # pedigree sample IDs are column 2 of the headerless, tab-delimited PED. keep.txt preserves the VCF's sample order
        awk -F'\\t' '{ sub(/\\r\$/, "") } \$2 != "" { print \$2 }' ${pedigree} | sort -u > ped_samples.txt
        bcftools query -l ${vcf} > vcf_samples.txt
        grep -Fxf ped_samples.txt vcf_samples.txt > keep.txt || true
        if [ ! -s keep.txt ]; then
            echo "ERROR: none of the pedigree samples are in ${vcf}" >&2
            exit 1
        fi
        if absent=\$(grep -vFxf vcf_samples.txt ped_samples.txt); then
            echo "WARNING: pedigree sample(s) absent from ${vcf}: \$(echo "\$absent" | paste -sd, -)" >&2
        fi

        read_region() {
            if [ "\$(wc -l < keep.txt)" -lt "\$(wc -l < vcf_samples.txt)" ]; then
                echo "Subsetting ${vcf} to \$(wc -l < keep.txt) of \$(wc -l < vcf_samples.txt) samples" >&2
                bcftools view ${region_args} ${subset_args} -S keep.txt -a -c 1 -Ou --no-version ${vcf} | \
                bcftools norm ${norm_args} -
            else
                bcftools norm ${norm_args} ${region_args} ${vcf}
            fi
        }

        read_region | \
        bcftools annotate \
            -x FORMAT/PL \
            -Ou \
            --no-version \
            - | \
        ${write_bcf}

        # -G drops all samples/genotypes - the annotation tools only need the sites
        bcftools view \
            -G \
            -Oz \
            --no-version \
            -o "${out_name}_sites.vcf.bgz" \
            "${out_name}_normalised.bcf"
        """
}
