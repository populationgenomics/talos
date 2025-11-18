"""
Add SpliceAi annotations from a reference HT
"""

from argparse import ArgumentParser

import loguru

import hail as hl
from hailtop.batch.job import BashJob

from cpg_utils import config, hail_batch

MISSING_FLOAT_LO = hl.float64(0.0)
MISSING_STRING = hl.str('missing')


def add_job(input_mt: str, output_mt: str, cohort_id: str) -> BashJob:
    """Adds a job to re-run this script."""

    job = hail_batch.get_batch().new_bash_job(f'Incorporate SpliceAi results for {cohort_id}')
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.command(f"""
    python -m talos.cpg_internal_scripts.cpgflow_jobs.AnnotateSpliceAiFromHt \\
        --input {input_mt} \\
        --output {output_mt}
    """)
    return job


def main(input_path: str, output_path: str):
    hail_batch.init_batch()

    mt = hl.read_matrix_table(input_path)

    ref_data_table = config.config_retrieve(['references', 'splice_ai_ht'])
    loguru.logger.info(f'Loading SpliceAi annotations from {input_path}')

    ht = hl.read_table(ref_data_table)

    # annotate the SpliceAi results into the main MatrixTable. Original schema:
    # 'splice_ai': struct {         noqa: ERA001
    #      delta_score: float32,    noqa: ERA001
    #      splice_consequence: str  noqa: ERA001
    #  }                            noqa: ERA001
    mt = mt.annotate_rows(splice_ai=ht[mt.row_key].splice_ai)

    # make sure these are real values, i.e. pad missings
    mt = mt.annotate_rows(
        splice_ai=mt.splice_ai.annotate(
            delta_score=hl.or_else(mt.splice_ai.delta_score, MISSING_FLOAT_LO),
            splice_consequence=hl.or_else(mt.splice_ai.splice_consequence, MISSING_STRING).replace(' ', '_'),
        ),
    )

    # print the resulting schema
    mt.describe()

    mt.write(output_path, overwrite=True)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', help='Path to the MatrixTable (input)')
    parser.add_argument('--output', help='Path to the MatrixTable (output)')
    args = parser.parse_args()

    main(args.input, args.output)
