process PrepareSvForSvAnnotate {
    container params.container

    input:
        tuple val(cohort), path(vcf), path(tbi), path(pedigree)

    // same VCF, minus any record with a non-symbolic ALT, with each complex variant's CPX_INTERVALS sorted
    output:
        tuple val(cohort), path("${cohort}_svannotate_ready.vcf.bgz"), path("${cohort}_svannotate_ready.vcf.bgz.tbi")

    script:
        // GATK SVAnnotate assumes CPX_INTERVALS is coordinate-sorted and aborts the whole run on the first
        // record where it is not - GATK-SV writes delINVdel as DEL,DEL,INV, which trips it immediately.
        // SVAnnotate also only works with symbolic alleles, so records with a non-symbolic ALT are removed here,
        // and each removed record is logged. Full reasoning lives in the script's docstring.
        // This runs ahead of AnnotateSvWithGatk rather than inside it because it needs the Talos container,
        // not the GATK one.
        //
        // The VCF is first cut down to the cohort's pedigree samples, dropping records (-c 1) no remaining sample
        // carries, so SVAnnotate and SVAFotate only see the cohort's SVs - same rules as NormaliseVcf. -I keeps the
        // INFO AC/AN as the source callset's: view would update those but not AF, leaving them inconsistent
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

        sv_input=${vcf}
        if [ "\$(wc -l < keep.txt)" -lt "\$(wc -l < vcf_samples.txt)" ]; then
            echo "Subsetting ${vcf} to \$(wc -l < keep.txt) of \$(wc -l < vcf_samples.txt) samples" >&2
            bcftools view -S keep.txt -I -c 1 -Ob --no-version -o ${cohort}_subset.bcf ${vcf}
            sv_input=${cohort}_subset.bcf
        fi

        python -m talos.scripts.prepare_sv_for_svannotate \
            --input \$sv_input \
            --output ${cohort}_svannotate_ready.vcf

        bgzip -c ${cohort}_svannotate_ready.vcf > ${cohort}_svannotate_ready.vcf.bgz
        tabix -p vcf ${cohort}_svannotate_ready.vcf.bgz
        """
}
