#!/usr/bin/env python3


"""
Simulated local run of Talos - Hail Query using a local runtime
Currently applies only to SV data as a minimal test case
"""

import os
from argparse import ArgumentParser
from datetime import datetime

from hailtop.batch.job import BashJob, Job

from cpg_utils.config import config_retrieve, output_path
from cpg_utils.hail_batch import authenticate_cloud_credentials_in_job, get_batch

from talos import (
    cpg_generate_pheno_ped,
    hail_filter_and_label,
    hail_filter_sv,
    hpo_panel_match,
    html_builder,
    query_panelapp,
    validate_categories,
)
from talos.static_values import get_logger


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


def pedigree_job(pedigree_in_gcp: str):
    """
    generates an extended pedigree from the Metamist content

    Returns:
        the pedigree file in the batch
    """
    new_job = get_batch().new_bash_job('Pedigree Generation')
    new_job.image(config_retrieve(['workflow', 'driver_image']))
    authenticate_cloud_credentials_in_job(new_job)
    new_job.command(f'{cpg_generate_pheno_ped.__file__} {config_retrieve(["workflow", "dataset"])} {new_job.output}')
    get_batch().write_output(new_job.output, pedigree_in_gcp)
    return new_job.output


def hpo_panel_job(ped_in_gcp: str, panel_file: str):
    """

    Args:
        ped_in_gcp ():
        panel_file ():

    Returns:

    """
    hpo_file = get_batch().read_input(config_retrieve(['workflow', 'obo_file']))
    hpo_job = get_batch().new_bash_job('Panel Matching')
    hpo_job.image(config_retrieve(['workflow', 'driver_image']))
    authenticate_cloud_credentials_in_job(hpo_job)
    hpo_job.command(f'python3 {hpo_panel_match.__file__} -i {ped_in_gcp} --hpo {hpo_file} --out {hpo_job.output}')
    get_batch().write_output(hpo_job.output, panel_file)
    return hpo_job.output


def panelapp_query_job(panel_file: str, panelapp_out: str):
    """

    Args:
        panel_file ():
        panelapp_out ():

    Returns:

    """
    panelapp_job = get_batch().new_bash_job('Panel Matching')
    panelapp_job.image(config_retrieve(['workflow', 'driver_image']))
    authenticate_cloud_credentials_in_job(panelapp_job)
    panelapp_job.command(
        f'python3 {query_panelapp.__file__} '
        f'--panels {panel_file} '
        f'--out_path {panelapp_job.output} '
        f'--dataset {config_retrieve(["workflow", "dataset"])}'
    )
    get_batch().write_output(panelapp_job.output, panelapp_out)
    return panelapp_job.output


def sort_out_sv(sv_path: str, panelapp: str, pedigree: str):
    sv_vcf_out = output_path('hail_SV_categorised.vcf.bgz', 'analysis')

    sv_job = get_batch().new_job('Local SV Filtering')
    set_job_resources(sv_job, cpu=2, memory='highmem')

    sv_job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

    mt_name = sv_path.split('/')[-1]

    sv_job.command('set -eux pipefail')

    # localise the input MT first (may take a while for chonky data)
    sv_job.command(f'gcloud --no-user-output-enabled storage cp -r {sv_path} .')
    # not sure if the other inputs need to be localised...
    sv_job.command(
        f'python3 {hail_filter_sv.__file__} --mt {mt_name} '
        f'--panelapp {panelapp} '
        f'--pedigree {pedigree} '
        f'--vcf_out {sv_job.output["vcf.bgz"]}',
    )
    get_batch().write_output(sv_job.output, sv_vcf_out.removesuffix('.vcf.bgz'))
    return sv_job.output


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
        cpu=8,
        memory='standard',
        storage=config_retrieve(['workflow', 'small_storage'], '500Gi'),
    )

    small_job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

    small_job.command('set -eux pipefail')

    clinvar_decisions = get_clinvar_table()
    clinvar_name = clinvar_decisions.split('/')[-1]

    # localise the clinvar decisions table
    small_job.command(
        f'cd $BATCH_TMPDIR && gcloud --no-user-output-enabled storage cp -r {clinvar_decisions} . && cd -'
    )

    # find, localise, and use the clinvar PM5 table
    pm5 = get_clinvar_table('clinvar_pm5')
    pm5_name = pm5.split('/')[-1]
    small_job.command(f'cd $BATCH_TMPDIR && gcloud --no-user-output-enabled storage cp -r {pm5} . && cd -')

    mt_name = mt_path.split('/')[-1]

    # now localise the input MT (may take a while for chonky data)
    small_job.command(f'cd ${{BATCH_TMPDIR}} && gcloud --no-user-output-enabled storage cp -r {mt_path} . && cd -')

    # not sure if the other inputs need to be localised...
    small_job.command(
        f'python3 {hail_filter_and_label.__file__} '
        f'--mt ${{BATCH_TMPDIR}}/{mt_name} '
        f'--panelapp {panelapp} '
        f'--pedigree {pedigree} '
        f'--vcf_out {small_job.output["vcf.bgz"]} '
        f'--checkpoint ${{BATCH_TMPDIR}}/checkpoint.mt '
        f'--clinvar ${{BATCH_TMPDIR}}/{clinvar_name} '
        f'--pm5 ${{BATCH_TMPDIR}}/{pm5_name} ',
    )
    get_batch().write_output(small_job.output, small_vcf_out.removesuffix('.vcf.bgz'))
    return small_job.output


def run_moi_tests(small_vcf, sv_vcf, json_results, panels, pedigree, party_panels):
    """
    run the MOI tests on the small and SV VCFs

    Args:
        small_vcf ():
        sv_vcf ():
        json_results ():
        panels ():
        pedigree ():
        party_panels ():
    """
    moi_job = get_batch().new_job('MOI Tests')
    set_job_resources(moi_job, memory='standard')

    moi_job.command(
        f'python3 {validate_categories.__file__} '
        f'--labelled_vcf {small_vcf} '
        f'--labelled_sv {sv_vcf} '
        f'--out_json {moi_job.output} '
        f'--panelapp {panels} '
        f'--pedigree {pedigree} '
        f'--participant_panels {party_panels}',
    )
    get_batch().write_output(moi_job.output, json_results)
    return moi_job.output


def make_html(moi_results, panels, html_output):
    """
    make an HTML

    Args:
        moi_results ():
        panels ():
        html_output ():
    """
    html_job = get_batch().new_job('Make HTML')
    set_job_resources(html_job, memory='standard')

    html_job.command(
        f'python3 {html_builder.__file__} '
        f'--results {moi_results} '
        f'--panelapp {panels} '
        f'--output {html_job.output} '
    )
    get_batch().write_output(html_job.output, html_output)


if __name__ == '__main__':
    get_logger(__name__).info(
        r"""Welcome To
 ███████████   █████████   █████          ███████     █████████
 █   ███   █  ███     ███   ███         ███     ███  ███     ███
     ███      ███     ███   ███        ███       ███ ███
     ███      ███████████   ███        ███       ███  █████████
     ███      ███     ███   ███        ███       ███         ███
     ███      ███     ███   ███      █  ███     ███  ███     ███
    █████    █████   █████ ███████████    ███████     █████████
        (local style)
        """,
    )

    parser = ArgumentParser()
    parser.add_argument('-i', help='Small variant data to analyse', required=True)
    parser.add_argument('-sv', help='SV data to analyse', required=False)
    args = parser.parse_args()

    # get pedigree
    ped_in_gcp = output_path('pedigree_extended.ped', 'analysis')
    ped_in_batch = pedigree_job(ped_in_gcp)

    # match HPOs
    panels_in_gcp = output_path('matched_panels.json', 'analysis')
    panel_file = hpo_panel_job(ped_in_batch, panels_in_gcp)

    # match panels
    panelapp_data = output_path('panelapp.json', 'analysis')
    panelapp_json = panelapp_query_job(panel_file, panelapp_data)

    # SVs in Hail
    sv_vcf = sort_out_sv(args.sv, panelapp_json, ped_in_batch)

    # small variants in Hail
    small_vcf = sort_out_smalls(args.i, panelapp_json, ped_in_batch)

    # results
    results_json = output_path('results.json', 'analysis')
    batch_results = run_moi_tests(small_vcf, sv_vcf, results_json, panelapp_json, ped_in_batch, panel_file)

    # HTML
    html_out = output_path('results.html', 'analysis')
    make_html(batch_results, panelapp_json, html_out)

    get_batch().run(wait=False)
