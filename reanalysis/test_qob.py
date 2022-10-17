#!/usr/bin/env python3

import os

from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, image_path, copy_common_env
from cpg_utils.git import prepare_git_job, get_organisation_name_from_current_directory, get_repo_name_from_current_directory, get_git_commit_ref_of_current_repository
import hailtop.batch as hb

QOB_SUB_SCRIPT = os.path.join(os.path.dirname(__file__), 'test_qob_sub.py')

def main():
    """
    main
    """

    service_backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch(
        name='AIP batch',
        backend=service_backend,
        cancel_after_n_failures=1,
        default_timeout=6000,
        default_memory='highmem',
    )

    j = batch.new_job().cpu(2).memory('highmem').storage('20G').image(image_path('hail'))

    prepare_git_job(
        job=j,
        organisation=get_organisation_name_from_current_directory(),
        repo_name=get_repo_name_from_current_directory(),
        commit=get_git_commit_ref_of_current_repository(),
    )
    copy_common_env(j)

    # might not be necessary, need to test.
    #j.command('sed -i s/gce/external/ /deploy-config/deploy-config.json')
    j.command('echo "Calling QoB test subscript."')
    j.command(f'python3 {QOB_SUB_SCRIPT}')

    batch.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=E1120


