process AnnotateVcfWithFreshClinvar {
    container params.container

    // apply this run's ClinVar decisions to one annotated shard. This is deliberately the last
    // annotation applied, and it happens in the per-run workflow so that a re-run always picks up
    // fresh ClinVar without invalidating the cached annotation half
    input:
        tuple val(cohort), path(vcf)
        path clinvar_zip

    output:
        tuple val(cohort), path("${vcf.simpleName}_clinvar.vcf.bgz")

    script:
        """
        set -euo pipefail

        echtvar anno \
            -e ${clinvar_zip} \
            ${vcf} \
            "${vcf.simpleName}_clinvar.vcf.bgz"
        """
}
