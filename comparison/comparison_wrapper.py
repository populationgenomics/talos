#!/usr/bin/env python3


"""
Entrypoint for the comparison process
"""

import logging
import os
import sys
from argparse import ArgumentParser

import hailtop.batch as hb

from cpg_utils.config import get_config
from cpg_utils.git import (
    get_git_commit_ref_of_current_repository,
    get_organisation_name_from_current_directory,
    get_repo_name_from_current_directory,
    prepare_git_job,
)
from cpg_utils.hail_batch import (
    authenticate_cloud_credentials_in_job,
    copy_common_env,
    image_path,
    output_path,
    remote_tmpdir,
)

# local script references
COMPARISON_SCRIPT = os.path.join(os.path.dirname(__file__), 'comparison.py')


def main(results_folder: str, seqr: str, mt: str):
    """
    main method, which runs the Talos comparison
    :param results_folder:
    :param seqr:
    :param mt:
    :return:
    """

    # set up a batch
    service_backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch(name='run Talos comparison', backend=service_backend)

    # create a new job
    comp_job = batch.new_job(name='Run Comparison')

    # set reasonable job resources
    comp_job.cpu(4).image(image_path('hail')).memory('standard').storage('50G')

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

    # need to localise the VCF + index
    run_vcf = os.path.join(results_folder, 'hail_categorised.vcf.bgz')
    vcf_in_batch = batch.read_input_group(**{'vcf.bgz': run_vcf, 'vcf.bgz.tbi': run_vcf + '.tbi'})
    ped_in_batch = batch.read_input(os.path.join(results_folder, 'latest_pedigree.fam'))

    results_command = (
        'pip install . && '
        f'python3 {COMPARISON_SCRIPT} '
        f'--results_folder {results_folder} '
        f'--pedigree {ped_in_batch} '
        f'--seqr {seqr} '
        f'--vcf {vcf_in_batch["vcf.bgz"]} '
        f'--mt {mt} '
        f'--output {output_path("comparison_result")} '
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
    parser = ArgumentParser()
    parser.add_argument('--results', help='results folder', required=True)
    parser.add_argument('--seqr', help='Flagged Seqr variants', required=True)
    parser.add_argument('--mt', help='Hail MT of annotated variants', required=True)
    args = parser.parse_args()
    main(results_folder=args.results, seqr=args.seqr, mt=args.mt)
