#!/usr/bin/env python3


"""
wrapper for reanalysis process

currently still reliant on a dataproc cluster for the Hail-VEP runtime environment

working to remove call out to Slivar,
compound-het calculations moving to Hail
"""


from typing import Any, Dict, Optional, Union

import json
import logging
import os
import sys

import click
from cloudpathlib import AnyPath
import hailtop.batch as hb

from analysis_runner import dataproc
from cpg_utils.hail import init_query_service, output_path

from query_panelapp import main as panelapp_main


DEFAULT_IMAGE = os.getenv('CPG_DRIVER_IMAGE')
assert DEFAULT_IMAGE

# static paths to write outputs
PANELAPP_JSON_OUT = output_path('panelapp_137_data.json')
HAIL_VCF_OUT = output_path('hail_categorised.vcf.bgz')
COMP_HET_JSON = output_path('hail_comp_het.json')
MT_OUT_PATH = output_path('hail_105_ac.mt')
CONFIG_OUT = output_path('config_used.json')
REHEADERED_OUT = output_path('hail_categories_reheadered.vcf.bgz')
MT_TMP = output_path('tmp_hail_table.mt', category='tmp')

# location of the Slivar Docker image
AR_REPO = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
BCFTOOLS_TAG = 'bcftools:1.10.2--h4f4756c_2'
BCFTOOLS_IMAGE = f'{AR_REPO}/{BCFTOOLS_TAG}'

# local script references
HAIL_SCRIPT = os.path.join(os.path.dirname(__file__), 'hail_filter_and_categorise.py')
RESULTS_SCRIPT = os.path.join(os.path.dirname(__file__), 'validate_categories.py')
DATAPROC_SETUP_SCRIPTS = [
    'gs://cpg-reference/hail_dataproc/install_common.sh',
    'gs://cpg-reference/vep/vep-GRCh38.sh',  # install & configure VEP 105
]


def read_json_dict_from_path(bucket_path: str) -> Dict[str, Any]:
    """
    take a path to a JSON file, read into an object
    :param bucket_path:
    """
    with open(AnyPath(bucket_path), encoding='utf-8') as handle:
        return json.load(handle)


def set_job_resources(job: Union[hb.batch.job.BashJob, hb.batch.job.Job]):
    """
    applied resources to the job
    :param job: apply resources to _this_ job
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

    script = (
        f'{HAIL_SCRIPT} '
        f'--mt {matrix} '
        f'--pap {PANELAPP_JSON_OUT} '
        f'--config {config} '
        f'--output {HAIL_VCF_OUT} '
        f'--mt_out {MT_OUT_PATH}'
        f'--mt_tmp {MT_TMP}'
    )
    hail_job = dataproc.hail_dataproc_job(
        batch=batch,
        worker_machine_type='n1-highmem-8',
        worker_boot_disk_size=200,
        secondary_worker_boot_disk_size=200,
        script=script,
        max_age='8h',
        init=DATAPROC_SETUP_SCRIPTS,
        job_name='run hail reanalysis stage',
        num_secondary_workers=10,
        num_workers=2,
        cluster_name='hail_reanalysis_stage',
    )
    set_job_resources(hail_job)

    return hail_job


def handle_reheader_job(
    batch: hb.Batch,
    local_vcf: str,
    config_dict: Dict[str, Any],
    prior_job: Optional[hb.batch.job.BashJob] = None,
) -> hb.batch.job.BashJob:
    """
    runs the bcftools re-header process
    :param batch:
    :param local_vcf:
    :param config_dict:
    :param prior_job:
    :return:
    """

    bcft_job = batch.new_job(name='bcftools_reheader_stage')
    set_job_resources(bcft_job)
    bcft_job.image(BCFTOOLS_IMAGE)

    if prior_job is not None:
        bcft_job.depends_on(prior_job)

    bcft_job.declare_resource_group(
        vcf={'vcf': '{root}.vcf.bgz', 'vcf.tbi': '{root}.vcf.bgz.tbi'}
    )

    # reheader the VCF using BCFtools and sed
    # replace the empty description with the full CSQ line from config
    desc = '##INFO=<ID=CSQ,Number=.,Type=String,Description="'

    # grotty string formatting to deliver the correct syntax to bcftools
    conf_csq = config_dict['variant_object'].get('csq_string').replace('|', r'\|')
    new_format = rf"Format: '{conf_csq}'"

    # only to reduce line length
    b_rh = 'bcftools reheader'

    bcft_job.command(
        'set -ex; '
        f'bcftools view -h {local_vcf} | sed \'s/'
        f'{desc}">/{desc}{new_format}">/\' > new_header; '
        f'{b_rh} -h new_header --threads 4 -o {bcft_job.vcf["vcf"]} {local_vcf}; '
        f'tabix {bcft_job.vcf["vcf"]}; '
    )
    return bcft_job


@click.command()
@click.option(
    '--matrix', 'matrix_path', help='variant matrix table to analyse', required=True
)
@click.option(
    '--config_json',
    help='dictionary of runtime settings',
    default='gs://cpg-acute-care-test/reanalysis/reanalysis_conf.json',
)
@click.option(
    '--panelapp_version',
    help='panelapp current comparison with this earlier version',
    required=False,
)
@click.option(
    '--panel_genes',
    help='location of a Gene list for use in analysis',
    required=False,
)
@click.option('--ped', 'ped_file', help='ped file for this analysis')
def main(
    matrix_path: str,
    config_json: str,
    panelapp_version: Optional[str],
    panel_genes: Optional[str],
    ped_file: str,
):
    """
    main method, which runs the full reanalysis process

    :param matrix_path:
    :param config_json:
    :param panelapp_version:
    :param panel_genes:
    :param ped_file:
    """

    logging.info('Starting the reanalysis batch')

    config_dict = read_json_dict_from_path(config_json)

    service_backend = hb.ServiceBackend(
        billing_project=os.getenv('HAIL_BILLING_PROJECT'),
        bucket=os.getenv('HAIL_BUCKET'),
    )

    # create a hail batch
    batch = hb.Batch(
        name='run_reanalysis', backend=service_backend, cancel_after_n_failures=1
    )

    # read ped and config files as a local batch resource
    # satisfy pylint until we use this in code
    print(ped_file)
    # ped_in_batch = batch.read_input(ped_file)

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
    # run hail categorisation #
    # ----------------------- #
    # permit continuity if the hail job isn't required
    prior_job = None

    # only run if the output VCF doesn't already exist
    if not AnyPath(HAIL_VCF_OUT).exists():
        prior_job = handle_hail_job(
            batch=batch,
            matrix=matrix_path,
            config=config_json,
        )

    # copy the Hail output file into the remaining batch jobs
    hail_output_in_batch = batch.read_input_group(
        **{'vcf': HAIL_VCF_OUT, 'vcf.tbi': HAIL_VCF_OUT + '.tbi'}
    )

    # --------------------------------- #
    # bcftools re-headering of hail VCF #
    # --------------------------------- #

    # this is no longer explicitly required...
    # it was required to run slivar: geneId, consequences, and transcript
    # if we can avoid using slivar for comp-hets, this isn't required
    # when extracting the consequences in python we can use the config string
    # this would mean the VCF is mostly useless when separated from the
    # config file... retain for now at least
    bcftools_job = handle_reheader_job(
        batch=batch,
        local_vcf=hail_output_in_batch['vcf'],
        config_dict=config_dict,
        prior_job=prior_job,
    )

    batch.write_output(bcftools_job.vcf, REHEADERED_OUT)

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
