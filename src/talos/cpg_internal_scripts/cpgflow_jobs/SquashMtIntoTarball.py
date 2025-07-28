from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def make_tarball_squash_job(
    cohort: targets.Cohort,
    input_mt: str,
    output_tar: Path,
    job_attrs: dict,
) -> 'BashJob':
    """
    Takes the MatrixTable, localises it, and squashes the deeply nested directory into a Tarball.
    """

    job = hail_batch.get_batch().new_job(
        name=f'CompressMtIntoTarball: {cohort.id}',
        attributes=job_attrs | {'tool': 'tar'},
    )
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.memory('highmem').cpu(4)
    job.storage(config.config_retrieve(['talos_prep', 'mt_storage'], '250Gi'))

    # copy the MT into the image, bundle it into a Tar-Ball
    job.command(f'gcloud --no-user-output-enabled storage cp --do-not-decompress -r {input_mt} $BATCH_TMPDIR')

    # rename the MT from cohort to dataset name - better flexibility downstream, probably harmonise this at some point
    job.command(f'mv $BATCH_TMPDIR/{cohort.id}.mt $BATCH_TMPDIR/{cohort.dataset.name}.mt')

    # once the data is copied - cd into the tmpdir, then tar it up
    job.command('cd $BATCH_TMPDIR')

    # no compression - the Hail objects are all internally gzipped so not much to gain there
    job.command(f'tar --remove-files -cf {job.output} {cohort.dataset.name}.mt')

    hail_batch.get_batch().write_output(job.output, output_tar)

    return job
