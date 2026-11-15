process EncodeClinvarEchtvar {
    container params.container
    publishDir params.processed_annotations, mode: 'copy'

    // encode the ClinvArbitration release VCF (all decisions, not just Pathogenic SNVs) into an
    // echtvar zip. The config renames the source INFO fields to the clinvar_* names the streaming
    // filter reads. Cohort-independent and month-stamped like the rest of ClinvArbitration's output,
    // so this runs once per month in the preparation workflow, not per Talos run
    input:
        path clinvar_vcf
        val timestamp

    output:
        path "clinvarbitration_${timestamp}.zip"

    script:
        """
        set -euo pipefail

        echtvar encode clinvarbitration_${timestamp}.zip /talos/echtvar/clinvar_config.json ${clinvar_vcf}
        """
}
