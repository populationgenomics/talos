from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def make_annotate_mito_job(
    cohort_id: str,
    mito_vcf: str,
    output: Path,
    job_attrs: dict,
) -> 'BashJob':
    """Annotate a Mito VCF with bcftools csq, and various mito annotations using echtvar."""

    batch_instance = hail_batch.get_batch('Talos')

    # localise the mito VCF, this is on chrM
    local_mito_vcf = batch_instance.read_input(mito_vcf)

    # get the GFF3 file required to generate consequences
    gff3_file = batch_instance.read_input(config.config_retrieve(['references', 'ensembl_gff3']))

    # get the fasta mito data was aligned/called on
    mito_fasta = config.config_retrieve(['references', 'mito_fasta'])
    fasta = batch_instance.read_input_group(fasta=mito_fasta, fai=f'{mito_fasta}.fai')['fasta']

    csq_job = batch_instance.new_bash_job(
        name=f'AnnotateMitoConsequenceWithBcftools: {cohort_id}',
        attributes=job_attrs | {'tool': 'bcftools'},
    )
    csq_job.image(config.config_retrieve(['workflow', 'driver_image']))
    csq_job.cpu(1).storage('10G')

    # then run the csq command:
    # -g is the GFF3 file
    # -B 10 is where to truncate long protein changes
    # -C 2 selects the vertebrate Mito codon lookup table
    # --local-csq indicates we want each variant annotated independently (haplotype unaware)
    # --unify-chr-names 'chr,-,chr' strips the `chr` prefix from the Fasta and VCF, to match the GFF file
    csq_job.command(
        f"""
        bcftools csq \\
            --force \\
            -f {fasta} \\
            --local-csq \\
            -C 2 \\
            --threads 4 \\
            -g {gff3_file} \\
            -B 10 \\
            --unify-chr-names 'chr,-,chr' \\
            -Oz -o bcftools_annotated.vcf.gz \\
            {local_mito_vcf}
        """,
    )

    # mito annotation sources formatted for echtvar
    mitimpact = batch_instance.read_input(config.config_retrieve(['references', 'mitimpact']))
    mitotip = batch_instance.read_input(config.config_retrieve(['references', 'mitotip']))
    napogee = batch_instance.read_input(config.config_retrieve(['references', 'napogee']))

    csq_job.command(
        f"""
        echtvar anno \
            -e ${napogee} \
            -e ${mitimpact} \
            -e ${mitotip} \
            bcftools_annotated.vcf.gz \
            {csq_job.output}
        """
    )

    batch_instance.write_output(csq_job.output, output)

    return csq_job
