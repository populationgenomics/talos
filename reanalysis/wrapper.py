#!/usr/bin/env python3


"""
wrapper for reanalysis process

currently still reliant on a dataproc cluster for the Hail-VEP runtime environment
"""


from typing import Any, Dict, Optional, Union
import json
import logging
import os
import sys
from shlex import quote
import click
from cloudpathlib import AnyPath
import hailtop.batch as hb

from analysis_runner import dataproc
from cpg_utils.hail import (
    init_query_service,
    output_path,
)

from reanalysis.query_panelapp import main as panelapp_main


# static paths to write outputs
PANELAPP_JSON_OUT = output_path('panelapp_137_data.json')
HAIL_VCF_OUT = output_path('hail_classified.vcf.bgz')
MT_OUT_PATH = output_path('hail_105_ac.mt')
CONFIG_OUT = output_path('config_used.json')


# location of the Slivar Docker image
AR_REPO = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
SLIVAR_TAG = 'slivar:v0.2.7'
BCFTOOLS_TAG = 'bcftools:1.10.2--h4f4756c_2'
SLIVAR_IMAGE = f'{AR_REPO}/{SLIVAR_TAG}'
BCFTOOLS_IMAGE = f'{AR_REPO}/{BCFTOOLS_TAG}'

# rubbish local references
HAIL_SCRIPT = os.path.join(os.path.dirname(__file__), 'hail_filter_and_classify.py')
RESULTS_SCRIPT = os.path.join(os.path.dirname(__file__), 'validate_classifications.py')


def read_json_dict_from_path(bucket_path: str) -> Dict[str, Any]:
    """
    take a GCP bucket path to a JSON file, read into an object
    this loop can read config files, or data
    :param bucket_path:
    :return:
    """
    with open(AnyPath(bucket_path), encoding='utf-8') as handle:
        return json.load(handle)


def set_job_resources(job: Union[hb.batch.job.BashJob, hb.batch.job.Job]):
    """
    applied resources to the job
    :param job:
    """
    job.cpu(2)
    job.memory('standard')
    job.storage('20G')


def handle_hail_job(batch: hb.Batch, matrix: str, config: str) -> hb.batch.job.Job:
    """
    sets up the hail-VEP105 dataproc environment
    :param batch:
    :param matrix:
    :param config:
    :return:
    """
    init_query_service()

    script = quote(
        f'{HAIL_SCRIPT} '
        f'--mt {matrix} '
        f'--pap {PANELAPP_JSON_OUT} '
        f'--config {config} '
        f'--output {HAIL_VCF_OUT} '
        f'--mt_out {MT_OUT_PATH}'
    )
    hail_job = dataproc.hail_dataproc_job(
        batch=batch,
        worker_machine_type='n1-highmem-8',
        worker_boot_disk_size=200,
        secondary_worker_boot_disk_size=200,
        script=script,
        max_age='8h',
        init=[
            'gs://cpg-reference/hail_dataproc/install_common.sh',
            'gs://cpg-reference/vep/vep-GRCh38.sh',  # install and configure VEP 105
        ],
        job_name='run hail reanalysis stage',
        num_secondary_workers=10,
        num_workers=2,
        cluster_name='hail_reanalysis_stage',
    )
    set_job_resources(hail_job)

    return hail_job


@click.command()
@click.option(
    '--matrix', 'matrix_path', help='variant matrix table to analyse', required=True
)
@click.option('--ped', 'ped_file', help='ped file for this analysis', required=True)
@click.option(
    '--config_json',
    help='dictionary of runtime settings',
    default='gs://cpg-acute-care-test/reanalysis/reanalysis_conf.json',
)
@click.option(
    '--panelapp_version',
    help='panelapp version for comparison with earlier version',
    required=False,
)
@click.option(
    '--panel_genes',
    help='location of a Gene list for use in analysis',
    required=False,
)
def main(
    matrix_path: str,
    config_json: str,
    panelapp_version: Optional[str],
    panel_genes: Optional[str],
):

    """
    main method, which runs the full reanalysis process

    :param matrix_path:
    :param config_json:
    :param panelapp_version:
    :param panel_genes:
    """

    logging.info('Starting the reanalysis batch')

    service_backend = hb.ServiceBackend(
        billing_project=os.getenv('HAIL_BILLING_PROJECT'),
        bucket=os.getenv('HAIL_BUCKET'),
    )

    # create a hail batch
    batch = hb.Batch(
        name='run_reanalysis', backend=service_backend, cancel_after_n_failures=1
    )

    # read config file as a local batch resource
    conf_in_batch = batch.read_input(config_json)

    # save a copy of the config file into the output location
    batch.write_output(conf_in_batch, CONFIG_OUT)

    # -------------------------------- #
    # query panelapp for panel details #
    # -------------------------------- #
    # no need to launch in a separate batch, minimal dependencies
    panelapp_main(
        panel_id='137',
        out_path=PANELAPP_JSON_OUT,
        previous_version=panelapp_version,
        gene_list=panel_genes,
    )

    # ----------------------- #
    # run hail classification #
    # ----------------------- #

    # only run if the output VCF doesn't already exist
    if not AnyPath(HAIL_VCF_OUT).exists():
        handle_hail_job(
            batch=batch,
            matrix=matrix_path,
            config=config_json,
        )

    # run the batch, and wait, so that the result metadata updates
    batch.run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%M-%d %H:%M:%S',
        stream=sys.stderr,
    )
    main()  # pylint: disable=E1120
