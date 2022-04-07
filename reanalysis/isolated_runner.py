#!/usr/bin/env python3


"""
test of the comp het process
"""


import logging
import os
import sys

import click
import hailtop.batch as hb

from analysis_runner.git import (
    prepare_git_job,
    get_repo_name_from_current_directory,
    get_git_commit_ref_of_current_repository,
)
from cpg_utils.hail import copy_common_env, output_path, remote_tmpdir

DEFAULT_IMAGE = os.getenv('CPG_DRIVER_IMAGE')
assert DEFAULT_IMAGE

# local script references
COMP_HET_SCRIPT = os.path.join(os.path.dirname(__file__), 'isolated_comphet_test.py')
COMP_HET_JSON = output_path('hail_comp_het.json')


@click.command()
@click.option('--matrix_path', help='variant matrix table to analyse', required=True)
def main(matrix_path: str):
    """
    main method, which runs the full reanalysis process

    :param matrix_path: annotated input matrix table
    """

    logging.info('Starting the reanalysis batch')

    service_backend = hb.ServiceBackend(
        billing_project=os.getenv('HAIL_BILLING_PROJECT'),
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch(
        name='run reanalysis (AIP)',
        backend=service_backend,
        cancel_after_n_failures=1,
    )

    labelling_job = batch.new_job(name='hail comp-het test')
    labelling_job.cpu(2).image(DEFAULT_IMAGE).memory('16Gi').storage('20G')
    prepare_git_job(
        job=labelling_job,
        repo_name=get_repo_name_from_current_directory(),
        commit=get_git_commit_ref_of_current_repository(),
    )
    labelling_command = (
        f'python3 {COMP_HET_SCRIPT} '
        f'--mt_input {matrix_path} '
        f'--out_json {COMP_HET_JSON} '
    )

    logging.info(f'PanelApp Command: {labelling_command}')
    labelling_job.command(labelling_command)
    copy_common_env(labelling_job)

    batch.run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )
    main()  # pylint: disable=E1120
