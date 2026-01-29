from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def make_condense_jobs(
    cohort_id: str,
    manifest_file: Path,
    manifest_dir: str,
    output: Path,
    tmp_dir: Path,
    job_attrs: dict,
) -> list['BashJob']:
    """read a manifest file, and generate a bash script to compose the VCF fragments into a single VCF file."""
    local_manifest = hail_batch.get_batch().read_input(manifest_file)

    # generate a bash script to do the composition
    job_1 = hail_batch.get_batch().new_bash_job(f'Create Compose Script: {cohort_id}', attributes=job_attrs)
    job_1.image(config.config_retrieve(['workflow', 'driver_image']))
    job_1.command(
        f"""
        python -m talos.cpg_internal_scripts.write_gcloud_compose_script \\
        --input {local_manifest} \\
        --vcf_dir {manifest_dir} \\
        --output {output!s} \\
        --script {job_1.output} \\
        --tmp {tmp_dir / 'compose_intermediates' / cohort_id!s}
        """,
    )

    job_2 = hail_batch.get_batch().new_bash_job(f'Run GCloud Compose: {cohort_id}', attributes=job_attrs)
    job_2.image(config.config_retrieve(['images', 'cpg_hail_gcloud']))
    job_2.command(f'bash {job_1.output}')

    return [job_1, job_2]
