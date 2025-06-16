from typing import TYPE_CHECKING

from cpg_utils import config, hail_batch, Path
from cpg_flow import targets


if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def make_annotation_transfer_job(
    dataset: targets.Dataset,
    annotations_ht: str,
    input_mt: str,
    output_mt: Path,
    job_attrs: dict,
) -> 'BashJob':
    """
    mixes the annotations HT with the variant MT, writing a final annotated MT.
    """

    job = hail_batch.get_batch().new_job(
        name=f'JumpAnnotationsFromHtToFinalMt: {dataset.name}',
        attributes=job_attrs | {'tool': 'hail'},
    )
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.cpu(16).memory('highmem').storage('250Gi')
    job.command(
        f"""
        python -m talos.annotation_scripts.TransferAnnotationsToMatrixTable \\
        --input_path {input_mt} \\
        --annotations {annotations_ht} \\
        --output {output_mt!s} \\
        --backend batch
        """,
    )

    return job
