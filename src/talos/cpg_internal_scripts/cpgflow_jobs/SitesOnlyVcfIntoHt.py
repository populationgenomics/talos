from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def make_vcf_to_ht_job(
    cohort_id: str,
    bcftools_vcf: str,
    output_ht: Path,
    tmp_dir: Path,
    job_attrs: dict,
) -> 'BashJob':
    """"""

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
        name=f'ProcessAnnotatedSitesOnlyVcfIntoHt: {cohort_id}',
        attributes=job_attrs | {'tool': 'hail'},
    )
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.cpu(2).storage('20Gi').memory('highmem')
    job.command(
        f"""
        python -m talos.annotation_scripts.ReformatAnnotatedVcfIntoHailTable \\
            --input {bcftools_vcf} \\
            --output {output_ht!s} \\
            --gene_bed {gene_roi} \\
            --am {alphamissense} \\
            --mane {mane_json} \\
            --checkpoint {tmp_dir!s}
        """,
    )

    return job
