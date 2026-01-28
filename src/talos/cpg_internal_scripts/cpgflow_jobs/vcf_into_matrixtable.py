from typing import TYPE_CHECKING


from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_mt_ingest_jobs(
    cohort_id: str,
    fragments: list[str],
    bcftools_template: str,
    checkpoint: str,
    outputs: dict[str, Path | str],
    job_attrs: dict[str, str],
):
    """One job per fragment, read each as MT and run some additional annotation/reformatting."""
    batch = hail_batch.get_batch()

    mane_json = batch.read_input(config.config_retrieve(['references', 'mane_json']))
    gene_roi = batch.read_input(config.config_retrieve(['references', 'ensembl_bed']))

    all_jobs: list['BashJob'] = []
    for part in fragments:
        vcf_path = bcftools_template.format(part=part)
        job = batch.new_bash_job(
            name=f'ConvertVcfToMt: {cohort_id}, {part}',
            attributes=job_attrs | {'tool': 'hail'},
        )
        output_mt = outputs['fragment_template'].format(part=part)

        # parallel checkpoint per fragment
        checkpoint_filled = checkpoint.format(part=part)

        job.image(config.config_retrieve(['workflow', 'driver_image']))
        job.cpu(2).storage('10G')
        job.command(f"""
        python -m talos.annotation_scripts \\
            --input {vcf_path} \\
            --output {output_mt} \\
            --gene_bed {gene_roi} \\
            --mane {mane_json} \\
            --checkpoint {checkpoint_filled}
        """)
        all_jobs.append(job)

    # create a final success job
    success = batch.new_bash_job(f'{cohort_id} bcftools success file')
    success.depends_on(*all_jobs)
    success.image(config.config_retrieve(['workflow', 'driver_image']))
    success.command(f'echo "Success!" > {success.output}')
    batch.write_output(success.output, outputs['success'])

    return [*all_jobs, success]
