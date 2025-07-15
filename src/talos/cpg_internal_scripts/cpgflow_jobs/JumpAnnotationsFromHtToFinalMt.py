from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def make_annotation_transfer_job(
    cohort_id: str,
    annotations_ht: str,
    input_mt: str,
    output_mt: Path,
    job_attrs: dict,
) -> 'BashJob':
    """
    mixes the annotations HT with the variant MT, writing a final annotated MT.
    """

    job = hail_batch.get_batch().new_job(
        name=f'JumpAnnotationsFromHtToFinalMt: {cohort_id}',
        attributes=job_attrs | {'tool': 'hail'},
    )
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    # using QOB, so no need for huge resources
    job.storage('10Gi')
    job.command(
        f"""
        python -m talos.annotation_scripts.TransferAnnotationsToMatrixTable \\
        --input {input_mt} \\
        --annotations {annotations_ht} \\
        --output {output_mt!s} \\
        --backend batch
        """,
    )

    return job
