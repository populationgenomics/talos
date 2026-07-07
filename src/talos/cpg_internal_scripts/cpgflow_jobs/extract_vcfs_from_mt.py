from random import randint
from typing import TYPE_CHECKING

from cpg_flow import targets, workflow
from cpg_utils import config, hail_batch, to_path

from talos.cpg_internal_scripts.cpg_flow_utils import query_for_latest_analysis

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def make_vcf_extraction_job(
    cohort: targets.Cohort,
    id_file: str,
    output: str,
    job_attrs: dict,
) -> 'BashJob':
    """Create a Hail Batch job to extract VCF from an AnnotateCohort MatrixTable."""

    # either get a mt from config, from metamist, or fail
    if not (
        input_mt := query_for_latest_analysis(
            dataset=workflow.get_multicohort().analysis_dataset.name,
            analysis_type='matrixtable',
            sequencing_type=config.config_retrieve(['workflow', 'sequencing_type']),
            long_read=config.config_retrieve(['workflow', 'long_read'], False),
            stage_name='AnnotateCohort',
        )
    ):
        raise ValueError(f'No MatrixTable found in Metamist for {cohort.id}')

    # write all SG IDs in this Cohort to a file
    with to_path(id_file).open('w') as f:
        for sg in cohort.get_sequencing_group_ids():
            f.write(f'{sg}\n')

    job = hail_batch.get_batch().new_job(f'ExtractDataFromMt: {cohort.id}', attributes=job_attrs)
    job.storage('10Gi')
    job.spot(False)
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    sgid_file_local = hail_batch.get_batch().read_input(id_file)
    job.command(f'sleep {randint(0, 1200)}')
    job.command(
        f"""
        python -m talos.cpg_internal_scripts.extract_fragmented_vcf_from_mt \\
            --input {input_mt} \\
            --sgs {sgid_file_local} \\
            --output {output!s}
        """,
    )

    return job
