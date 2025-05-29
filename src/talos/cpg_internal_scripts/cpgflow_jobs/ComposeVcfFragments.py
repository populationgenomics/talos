from typing import TYPE_CHECKING

from cpg_utils import config, hail_batch, Path
from cpg_flow import targets


if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def make_condense_jobs(
    dataset: targets.Dataset,
    manifest_file: str,
    output: Path,
    job_attrs: dict,
) -> list['BashJob']:
    local_manifest = hail_batch.get_batch().read_input(manifest_file)

    # generate a bash script to do the composition
    job_1 = hail_batch.get_batch().new_bash_job(f'Create Compose Script: {dataset.name}', attributes=job_attrs)
    job_1.image(config.config_retrieve(['workflow', 'driver_image']))
    job_1.command(
        f"""
        python -m talos.cpg_internal_scripts.write_gcloud_compose_script \\
        --input {local_manifest} \\
        --output {output!s} \\
        --script {job_1.output} \\
        --tmp
        """,
    )

    job_2 = hail_batch.get_batch().new_bash_job(f'Run GCloud Compose: {dataset.name}', attributes=job_attrs)
    job_2.image(config.config_retrieve(['images', 'cpg_hail_gcloud']))
    job_2.command(f'bash {job_1.output}')

    return [job_1, job_2]
