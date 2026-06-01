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
) -> 'BashJob':
    """read a manifest file, and generate a bash script to compose the VCF fragments into a single VCF file."""
    local_manifest = hail_batch.get_batch().read_input(manifest_file)

    # generate a bash script to do the composition
    job = hail_batch.get_batch().new_bash_job(f'Create & Run Compose Script: {cohort_id}', attributes=job_attrs)
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.command(
        f"""
        python -m talos.cpg_internal_scripts.write_gcloud_compose_script \\
        --input {local_manifest} \\
        --vcf_dir {manifest_dir} \\
        --output {output!s} \\
        --script condense_script.sh \\
        --tmp {tmp_dir / 'compose_intermediates' / cohort_id!s}

        bash condense_script.sh
        """,
    )

    return job
