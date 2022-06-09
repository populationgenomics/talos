#!/usr/bin/env python3


"""
Entrypoint for the comparison process
initially just print results
"""


import logging
import os
import sys

import click
import hailtop.batch as hb

from cpg_utils.git import (
    prepare_git_job,
    get_git_commit_ref_of_current_repository,
    get_organisation_name_from_current_directory,
    get_repo_name_from_current_directory,
)
from cpg_utils.hail_batch import (
    authenticate_cloud_credentials_in_job,
    copy_common_env,
    remote_tmpdir,
)
from cpg_utils.config import get_config


# local script references
COMPARISON_SCRIPT = os.path.join(os.path.dirname(__file__), 'comparison.py')


@click.command()
@click.option('--results', help='AIP result JSON')
@click.option('--seqr', help='Seqr flagged variants export')
@click.option('--ped', help='plink file for the cohort')
def main(results: str, seqr: str, ped: str):
    """
    main method, which runs the AIP comparison
    :param results:
    :param seqr:
    :param ped:
    :return:
    """

    # set up a batch
    service_backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch(
        name='run AIP comparison',
        backend=service_backend,
        cancel_after_n_failures=1,
        default_timeout=1000,
    )

    # create a new job
    comp_job = batch.new_job(name='Run Comparison')

    image_default = get_config()['workflow']['driver_image']

    # set reasonable job resources
    comp_job.cpu(2).image(image_default).memory('standard').storage('50G')

    # run gcloud authentication
    authenticate_cloud_credentials_in_job(comp_job)

    # copy in Env Variables from current config
    copy_common_env(comp_job)

    # copy the relevant scripts into a Driver container instance
    prepare_git_job(
        job=comp_job,
        organisation=get_organisation_name_from_current_directory(),
        repo_name=get_repo_name_from_current_directory(),
        commit=get_git_commit_ref_of_current_repository(),
    )

    results_command = (
        'pip install peddy==0.4.8 && '
        f'PYTHONPATH=$(pwd) python3 {COMPARISON_SCRIPT} '
        f'--results {results} '
        f'--seqr {seqr} '
        f'--ped {ped} '
    )
    logging.info(f'Results command: {results_command}')
    comp_job.command(results_command)

    batch.run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )
    main()  # pylint: disable=E1120
