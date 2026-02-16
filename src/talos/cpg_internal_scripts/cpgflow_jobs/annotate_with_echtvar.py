from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def make_echtvar_job(
    cohort_id: str,
    fragments: dict[str, Path],
    am_zip: str,
    outputs: dict[str, Path],
    use_avi: bool,
    job_attrs: dict,
) -> list['BashJob']:
    """Create a job to annotate gnomAD frequencies and AlphaMissense using the Echtvar tool."""

    batch = hail_batch.get_batch()

    gnomad_annotations = batch.read_input(config.config_retrieve(['references', 'echtvar_gnomad']))

    avi_command = ''
    if use_avi:
        avi_annotations = batch.read_input(config.config_retrieve(['references', 'avi_zip']))
        avi_command = f'-e {avi_annotations}'

    # collect all jobs
    all_jobs: list['BashJob'] = []  # noqa: UP037

    for part, vcf_path in fragments.items():
        local_vcf = batch.read_input_group(
            vcf=vcf_path,
            vcf_idx=f'{vcf_path}.tbi',
        ).vcf

        am_local = batch.read_input(am_zip)

        job = batch.new_job(
            name=f'AnnotateWithEchtvar: {cohort_id}, {part}',
            attributes=job_attrs | {'tool': 'echtvar'},
        )
        job.image(config.config_retrieve(['images', 'echtvar']))
        job.command(f"""
        echtvar anno \\
            -e {gnomad_annotations} {avi_command} \\
            -e {am_local} \\
            -i "gnomad_AF_joint < 0.05" \\
            {local_vcf} \\
            {job.output}
        """)

        job.storage('50Gi').memory('highmem').cpu(4)
        batch.write_output(job.output, str(outputs['fragment_template']).format(part=part))
        all_jobs.append(job)

    # create a final success job
    success = batch.new_bash_job(f'{cohort_id} echtvar success file')
    success.depends_on(*all_jobs)
    success.image(config.config_retrieve(['workflow', 'driver_image']))
    success.command(f'echo "Success!" > {success.output}')
    batch.write_output(success.output, outputs['success'])

    return [*all_jobs, success]
