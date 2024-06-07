#!/usr/bin/env python3


"""
Simulated local run of AIP - Hail Query using a local runtime
Currently applies only to SV data as a minimal test case
"""

import os
from argparse import ArgumentParser
from datetime import datetime

from hailtop.batch.job import BashJob, Job

from cpg_utils.config import config_retrieve, output_path
from cpg_utils.hail_batch import authenticate_cloud_credentials_in_job, get_batch

from reanalysis import hail_filter_and_label, hail_filter_sv
from reanalysis.static_values import get_logger


def get_clinvar_table(key: str = 'clinvar_decisions') -> str:
    """
    try and identify the clinvar table to use
    - try the config specified path
    - fall back to storage:common default path
    - failing that, stick to standard annotations

    Args
        key (str): the key to look for in the config

    Returns:
        a path to a clinvar table, or None
    """

    if (clinvar_table := config_retrieve(['workflow', key], None)) is not None:
        get_logger().info(f'Using clinvar table {clinvar_table} from config')
        return clinvar_table

    get_logger().info(f'No forced {key} table available, trying default')

    clinvar_table = os.path.join(
        config_retrieve(['storage', 'common', 'analysis']),
        'aip_clinvar',
        datetime.now().strftime('%y-%m'),
        f'{key}.ht',
    )

    get_logger().info(f'Using clinvar table {clinvar_table}')
    return clinvar_table


def set_job_resources(job: BashJob, cpu: int = 4, memory: str = 'highmem', storage: str = '10Gi'):
    """
    apply resources to the job

    Args:
        job (Job): the job to set resources on
        cpu (int): number of CPUs to use
        memory (str): lowmem/standard/highmem
        storage (str): storage setting to use
    """
    # apply all settings to this job
    job.cpu(cpu).image(config_retrieve(['workflow', 'driver_image'])).memory(memory).storage(storage)
    authenticate_cloud_credentials_in_job(job)


def sort_out_sv(sv_path: str, panelapp: str, pedigree: str):
    sv_vcf_out = output_path('hail_SV_categorised.vcf.bgz', 'analysis')

    sv_job = get_batch().new_job('Local SV Filtering')
    set_job_resources(sv_job, cpu=2, memory='highmem')

    sv_job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

    mt_name = sv_path.split('/')[-1]

    sv_job.command('set -eux pipefail')

    # localise the input MT first (may take a while for chonky data)
    sv_job.command(f'gcloud --verbosity=none storage cp -r {sv_path} .')
    # not sure if the other inputs need to be localised...
    sv_job.command(
        f'python3 {hail_filter_sv.__file__} --mt {mt_name} '
        f'--panelapp {get_batch().read_input(panelapp)} '
        f'--pedigree {get_batch().read_input(pedigree)} '
        f'--vcf_out {sv_job.output["vcf.bgz"]}',
    )
    get_batch().write_output(sv_job.output, sv_vcf_out.removesuffix('.vcf.bgz'))


def sort_out_smalls(mt_path: str, panelapp: str, pedigree: str):
    """
    run the small variant classification and filtering on a local MT
    this needs the relevant ClinVar tables to be copied in

    Args:
        mt_path ():
        panelapp ():
        pedigree ():
    """
    small_vcf_out = output_path('hail_small_categorised.vcf.bgz', 'analysis')

    small_job = get_batch().new_job('Local Small Variant Filtering')
    set_job_resources(
        small_job,
        cpu=16,
        memory='standard',
        storage=config_retrieve(['workflow', 'small_storage'], '500Gi'),
    )

    small_job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

    small_job.command('set -eux pipefail')

    clinvar_decisions = get_clinvar_table()
    clinvar_name = clinvar_decisions.split('/')[-1]

    # localise the clinvar decisions table
    small_job.command(f'cd $BATCH_TMPDIR && gcloud storage cp -r {clinvar_decisions} . && cd -')

    # find, localise, and use the clinvar PM5 table
    pm5 = get_clinvar_table('clinvar_pm5')
    pm5_name = pm5.split('/')[-1]
    small_job.command(f'cd $BATCH_TMPDIR && gcloud storage cp -r {pm5} . && cd -')

    mt_name = mt_path.split('/')[-1]

    # now localise the input MT (may take a while for chonky data)
    small_job.command(f'cd ${{BATCH_TMPDIR}} && gcloud storage cp -r {mt_path} . && cd -')

    # not sure if the other inputs need to be localised...
    small_job.command(
        f'python3 {hail_filter_and_label.__file__} '
        f'--mt ${{BATCH_TMPDIR}}/{mt_name} '
        f'--panelapp {get_batch().read_input(panelapp)} '
        f'--pedigree {get_batch().read_input(pedigree)} '
        f'--vcf_out {small_job.output["vcf.bgz"]} '
        f'--checkpoint ${{BATCH_TMPDIR}}/checkpoint.mt '
        f'--clinvar ${{BATCH_TMPDIR}}/{clinvar_name} '
        f'--pm5 ${{BATCH_TMPDIR}}/{pm5_name} ',
    )
    get_batch().write_output(small_job.output, small_vcf_out.removesuffix('.vcf.bgz'))


if __name__ == '__main__':
    get_logger(__name__).info(
        r"""Welcome To The
          ___  _____ ______
         / _ \|_   _|| ___ \
        / /_\ \ | |  | |_/ /
        |  _  | | |  |  __/
        | | | |_| |_ | |
        \_| |_/\___/ \_|
        (local style)
        """,
    )

    parser = ArgumentParser()
    parser.add_argument('-i', help='Small variant data to analyse', required=True)
    parser.add_argument('-sv', help='SV data to analyse', required=False)
    parser.add_argument('--pedigree', help='PED file for the cohort', required=True)
    parser.add_argument('--panels', help='PanelApp Results!')
    args = parser.parse_args()

    sort_out_sv(args.sv, args.panels, args.pedigree)
    sort_out_smalls(args.i, args.panels, args.pedigree)

    get_batch().run(wait=False)
