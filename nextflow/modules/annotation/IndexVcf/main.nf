process IndexVcf {
    container params.container

    // single-VCF inputs may arrive without a .tbi; region scatter needs one. Only the index is
    // emitted (the VCF is joined back in the workflow) so the input is never re-uploaded as output
    input:
        tuple val(cohort), path(vcf)

    output:
        tuple val(cohort), path("${vcf}.tbi")

    script:
        """
        set -euo pipefail

        tabix ${vcf}
        """
}
