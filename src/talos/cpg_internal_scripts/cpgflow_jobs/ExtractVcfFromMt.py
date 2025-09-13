from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_utils import Path, config, hail_batch

from talos.cpg_internal_scripts.cpg_flow_utils import query_for_latest_analysis

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def make_vcf_extraction_job(
    cohort: targets.Cohort,
    output_mt: Path,
    output_sitesonly: Path,
    job_attrs: dict,
) -> 'BashJob':
    """Create a Hail Batch job to extract VCF from a dataset MatrixTable."""

    # either get a mt from config, from metamist, or fail
    if not (input_mt := config.config_retrieve(['workflow', 'starting_mt'], None)):  # noqa: SIM102
        # grab one from metamist instead, parameterised by sequencing type
        if not (
            input_mt := query_for_latest_analysis(
                dataset=cohort.dataset.name,
                analysis_type='matrixtable',
                sequencing_type=config.config_retrieve(['workflow', 'sequencing_type']),
                long_read=config.config_retrieve(['workflow', 'long_read'], False),
            )
        ):
            raise ValueError(f'No MatrixTable found in Metamist for {cohort.id}')

    # get the BED file - does not need to be localised
    bed = config.config_retrieve(['references', 'ensembl_merged_bed'])

    job = hail_batch.get_batch().new_job(f'ExtractDataFromDatasetMt: {cohort.id}', attributes=job_attrs)
    job.storage('10Gi')
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.command(
        f"""
        python -m talos.cpg_internal_scripts.extract_fragmented_vcf_from_mt \\
            --input {input_mt} \\
            --output_mt {output_mt!s} \\
            --output_sites_only {output_sitesonly!s} \\
            --bed {bed}
        """,
    )

    return job
