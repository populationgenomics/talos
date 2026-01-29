from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def make_bcftools_anno_jobs(
    cohort_id: str,
    fragments: list[str],
    echtvar_template: str,
    outputs: dict[str, Path | str],
    job_attrs: dict,
) -> list['BashJob']:
    """Use BCFtools to annotate each VCF fragment using CSQ"""
    batch = hail_batch.get_batch()

    gff3_file = hail_batch.get_batch().read_input(config.config_retrieve(['references', 'ensembl_gff3']))
    fasta = hail_batch.get_batch().read_input(config.config_retrieve(['references', 'ref_fasta']))

    all_jobs: list['BashJob'] = []  # noqa: UP037
    for part in fragments:
        vcf_path = echtvar_template.format(part=part)
        local_vcf = batch.read_input_group(
            vcf=vcf_path,
            vcf_idx=f'{vcf_path}.tbi',
        ).vcf

        job = hail_batch.get_batch().new_bash_job(
            name=f'AnnotateConsequenceWithBcftools: {cohort_id}, {part}',
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
            bcftools index -t {local_vcf}
            bcftools csq --force -f {fasta} \\
                --local-csq \\
                --threads 4 \\
                -g {gff3_file} \\
                -B 10 \\
                --unify-chr-names 'chr,-,chr' \\
                -Oz -o {job.output['vcf.bgz']} \\
                -W=tbi \\
                {local_vcf}
            """,
        )

        hail_batch.get_batch().write_output(
            job.output,
            str(outputs['fragment_template']).format(part=part).removesuffix('.vcf.bgz'),
        )
        all_jobs.append(job)

    # create a final success job
    success = batch.new_bash_job(f'{cohort_id} bcftools success file')
    success.depends_on(*all_jobs)
    success.image(config.config_retrieve(['workflow', 'driver_image']))
    success.command(f'echo "Success!" > {success.output}')
    batch.write_output(success.output, outputs['success'])

    return [*all_jobs, success]
