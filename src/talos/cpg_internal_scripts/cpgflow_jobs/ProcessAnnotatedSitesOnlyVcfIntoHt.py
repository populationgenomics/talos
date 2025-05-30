from typing import TYPE_CHECKING

from cpg_utils import config, hail_batch, Path
from cpg_flow import targets


if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def make_vcf_to_ht_job(
        dataset: targets.Dataset,
        bcftools_vcf: str,
        output_ht: Path,
        tmp_dir: Path,
        job_attrs: dict,
):
    """

    Args:
        dataset ():
        bcftools_vcf ():
        output_ht ():
        tmp_dir ():
        job_attrs ():

    Returns:

    """

    # pull the alphamissense TarBall location from config, and localise it
    alphamissense = config.reference_path('alphamissense/ht')

    # mane version for gene details
    mane_version = config.config_retrieve(['workflow', 'mane_version'], '1.4')
    mane_json = hail_batch.get_batch().read_input(config.reference_path(f'mane_{mane_version}/json'))

    # ensembl version used to generate region of interest
    ensembl_version = config.config_retrieve(['workflow', 'ensembl_version'], 113)
    # local file parsed into a dict
    gene_roi = hail_batch.get_batch().read_input(config.reference_path(f'ensembl_{ensembl_version}/bed'))

    job = hail_batch.get_batch().new_job(
        name=f'ProcessAnnotatedSitesOnlyVcfIntoHt: {dataset.name}',
        attributes=job_attrs | {'tool': 'hail'},
    )
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.cpu(4).storage('20Gi').memory('highmem')
    job.command(
        f'convert_annotated_vcf_to_ht '
        f'--input {bcftools_vcf} '
        f'--output {output_ht!s} '
        f'--gene_bed {gene_roi} '
        f'--am {alphamissense} '
        f'--mane {mane_json} '
        f'--checkpoint_dir {tmp_dir!s} ',
    )

    return job
