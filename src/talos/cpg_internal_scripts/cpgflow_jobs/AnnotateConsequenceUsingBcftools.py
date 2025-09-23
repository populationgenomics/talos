from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def make_bcftools_anno_job(
    cohort_id: str,
    gnomad_vcf: str,
    output: Path,
    job_attrs: dict,
) -> 'BashJob':
    """"""

    gnomad_annotated_vcf = hail_batch.get_batch().read_input(gnomad_vcf)

    # get the GFF3 file required to generate consequences
    gff3_file = hail_batch.get_batch().read_input(config.config_retrieve(['references', 'ensembl_gff3']))

    # get the fasta
    fasta = hail_batch.get_batch().read_input(config.config_retrieve(['references', 'ref_fasta']))

    job = hail_batch.get_batch().new_bash_job(
        name=f'AnnotateConsequenceWithBcftools: {cohort_id}',
        attributes=job_attrs | {'tool': 'bcftools'},
    )
    job.image(config.config_retrieve(['images', 'bcftools']))
    job.cpu(4).storage('20G')

    job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

    # the echtvar image doesn't have a tool to index, so first add that
    # then run the csq command:
    # -g is the GFF3 file
    # -B 10 is where to truncate long protein changes
    # --local-csq indicates we want each variant annotated independently (haplotype unaware)
    job.command(
        f"""
        bcftools index -t {gnomad_annotated_vcf}
        bcftools csq --force -f {fasta} \\
            --local-csq \\
            --threads 4 \\
            -g {gff3_file} \\
            -B 10 \\
            -Oz -o {job.output['vcf.bgz']} \\
            -W=tbi \\
            {gnomad_annotated_vcf}
        """,
    )

    hail_batch.get_batch().write_output(job.output, str(output).removesuffix('.vcf.bgz'))

    return job
