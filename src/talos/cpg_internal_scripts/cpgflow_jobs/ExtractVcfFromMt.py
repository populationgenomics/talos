from typing import TYPE_CHECKING

from cpg_utils import config, hail_batch, Path
from cpg_flow import targets

from talos.cpg_internal_scripts.cpg_flow_utils import query_for_latest_analysis


if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def make_vcf_extraction_job(
    dataset: targets.Dataset,
    output_mt: Path,
    output_sitesonly: Path,
    job_attrs: dict,
) -> 'BashJob':
    """
    Create a Hail Batch job to extract VCF from a dataset MatrixTable.
    """

    # either get a mt from metamist, or take one from config
    if not (
        input_mt := query_for_latest_analysis(
            dataset=dataset.name,
            analysis_type='matrixtable',
            object_suffix='.mt',
        )
    ):
        raise ValueError(f'No MatrixTable found in Metamist for {dataset.name}')

    # get the BED file - does not need to be localised
    ensembl_version = config.config_retrieve(['workflow', 'ensembl_version'], 113)
    bed = config.reference_path(f'ensembl_{ensembl_version}/merged_bed')
    job = hail_batch.get_batch().new_job(f'ExtractDataFromDatasetMt: {dataset.name}', attributes=job_attrs)
    job.storage('10Gi')
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.command(
        f"""
        python -m talos.cpg_internal_scripts.extract_fragmented_vcf_from_mt \
        --input {input_mt} \
        --output_mt {output_mt!s} \
        --output_sites_only {output_sitesonly!s} \
        --bed {bed}
        """,
    )

    return job
