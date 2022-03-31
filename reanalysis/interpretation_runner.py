#!/usr/bin/env python3


"""
Entrypoint for the interpretation pipeline process, runs the end-to-end
pipeline stages either directly or via Hail Batch(es)
 - Data extraction from PanelApp
 - Filtering and Annotation of variant data
 - Re-headering of resultant VCF

Steps are run only where the specified output does not exist
i.e. the full path to the output file is crucial, and forcing steps to
re-run currently requires the deletion of previous outputs

compound-het calculations moved to Hail, removed requirement for Slivar stage
"""


from typing import Any, Dict, Optional, Union

import json
import logging
import os
import sys

import click
from cloudpathlib import AnyPath
import hailtop.batch as hb

from analysis_runner.git import (
    prepare_git_job,
    get_repo_name_from_current_directory,
    get_git_commit_ref_of_current_repository,
)
from cpg_utils.hail import copy_common_env, image_path, output_path, remote_tmpdir


# static paths to write outputs
PANELAPP_JSON_OUT = output_path('panelapp_137_data.json')
HAIL_VCF_OUT = output_path('hail_categorised.vcf.bgz')
COMP_HET_JSON = output_path('hail_comp_het.json')
REHEADERED_OUT = output_path('hail_categories_reheadered.vcf.bgz')
REHEADERED_PREFIX = output_path('hail_categories_reheadered')
MT_TMP = output_path('tmp_hail_table.mt', category='tmp')
RESULTS_JSON = output_path('summary_results.json')

# location of the CPG BCFTools image
BCFTOOLS_IMAGE = image_path('bcftools:1.10.2--h4f4756c_2')
DEFAULT_IMAGE = os.getenv('CPG_DRIVER_IMAGE')
assert DEFAULT_IMAGE

# local script references
HAIL_FILTER = os.path.join(os.path.dirname(__file__), 'hail_filter_and_label.py')
QUERY_PANELAPP = os.path.join(os.path.dirname(__file__), 'query_panelapp.py')
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


def set_job_resources(
    job: Union[hb.batch.job.BashJob, hb.batch.job.Job],
    git=False,
    image: Optional[str] = None,
    prior_job: Optional[hb.batch.job.Job] = None,
):
    """
    applied resources to the job
    :param job:
    :param git:
    :param image:
    :param prior_job:
    """
    job.cpu(2)
    job.memory('standard')
    job.storage('20G')
    if prior_job is not None:
        job.depends_on(prior_job)

    job.image(image or DEFAULT_IMAGE)
    if git:
        # copy the relevant scripts into a Driver container instance
        prepare_git_job(
            job=job,
            repo_name=get_repo_name_from_current_directory(),
            commit=get_git_commit_ref_of_current_repository(),
        )


def handle_panelapp_job(
    batch: hb.Batch,
    gene_list: Optional[str],
    prev_version: Optional[str],
    prior_job: Optional[hb.batch.job.Job] = None,
) -> hb.batch.job.Job:
    """

    :param batch:
    :param gene_list:
    :param prev_version:
    :param prior_job:
    """
    panelapp_job = batch.new_job(name='query panelapp')
    set_job_resources(panelapp_job, git=True, prior_job=prior_job)
    panelapp_command = (
        f'python3 {QUERY_PANELAPP} --panel_id 137 --out_path {PANELAPP_JSON_OUT} '
    )
    if gene_list is not None:
        panelapp_command += f'--gene_list {gene_list} '
    elif prev_version is not None:
        panelapp_command += f'--previous_version {prev_version} '

    logging.info(f'PanelApp Command: {panelapp_command}')
    panelapp_job.command(panelapp_command)
    return panelapp_job


def handle_hail_filtering(
    batch: hb.Batch,
    matrix_path: str,
    config: str,
    prior_job: Optional[hb.batch.job.Job] = None,
) -> hb.batch.job.BashJob:
    """
    hail-query backend version of the filtering implementation
    use the init query service instead of running inside dataproc

    :param batch:
    :param matrix_path: path to annotated matrix table
    :param config:
    :param prior_job:
    :return:
    """

    labelling_job = batch.new_job(name='hail filtering')
    set_job_resources(labelling_job, git=True, prior_job=prior_job)
    labelling_command = (
        f'python3 {HAIL_FILTER} '
        f'--mt_input {matrix_path} '
        f'--panelapp_path {PANELAPP_JSON_OUT} '
        f'--config_path {config} '
        f'--out_vcf {HAIL_VCF_OUT} '
        f'--mt_tmp {MT_TMP}'
    )

    logging.info(f'PanelApp Command: {labelling_command}')
    labelling_job.command(labelling_command)
    copy_common_env(labelling_job)
    return labelling_job


def handle_reheader_job(
    batch: hb.Batch,
    local_vcf: str,
    config_dict: Dict[str, Any],
    prior_job: Optional[hb.batch.job.Job] = None,
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
    set_job_resources(bcft_job, image=BCFTOOLS_IMAGE, prior_job=prior_job)

    bcft_job.declare_resource_group(
        vcf={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'}
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
        f'{b_rh} -h new_header --threads 4 -o {bcft_job.vcf["vcf.bgz"]} {local_vcf}; '
        f'tabix {bcft_job.vcf["vcf.bgz"]}; '
    )
    return bcft_job


def handle_results_job(
    batch: hb.Batch,
    config: str,
    comp_het: str,
    reheadered_vcf: str,
    pedigree: str,
    prior_job: Optional[hb.batch.job.Job] = None,
) -> hb.batch.job.Job:
    """

    :param batch:
    :param config:
    :param comp_het:
    :param reheadered_vcf:
    :param pedigree:
    :param prior_job:
    :return:
    """
    results_job = batch.new_job(name='finalise_results')
    set_job_resources(results_job, git=True, prior_job=prior_job)
    results_command = (
        'pip install cyvcf2==0.30.14 && '
        f'PYTHONPATH=$(pwd) python3 {RESULTS_SCRIPT} '
        f'--config_path {config} '
        f'--comp_het {comp_het} '
        f'--class_vcf {reheadered_vcf} '
        f'--panelapp {PANELAPP_JSON_OUT} '
        f'--pedigree {pedigree} '
        f'--out_json {RESULTS_JSON} '
    )
    logging.info(f'Results command: {results_command}')
    results_job.command(results_command)
    return results_job


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
@click.option('--pedigree', help='location of a PED file')
def main(
    matrix_path: str,
    config_json: str,
    panelapp_version: Optional[str],
    panel_genes: Optional[str],
    pedigree: str,
):
    """
    main method, which runs the full reanalysis process

    :param matrix_path: annotated input matrix table
    :param config_json:
    :param panelapp_version:
    :param panel_genes:
    :param pedigree:
    """

    logging.info('Starting the reanalysis batch')

    config_dict = read_json_dict_from_path(config_json)

    service_backend = hb.ServiceBackend(
        billing_project=os.getenv('HAIL_BILLING_PROJECT'),
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch(
        name='run reanalysis (AIP)',
        backend=service_backend,
        cancel_after_n_failures=1,
    )

    # -------------------------------- #
    # query panelapp for panel details #
    # -------------------------------- #
    # no need to launch in a separate batch, minimal dependencies
    prior_job = handle_panelapp_job(
        batch=batch, gene_list=panel_genes, prev_version=panelapp_version
    )

    # ----------------------- #
    # run hail categorisation #
    # ----------------------- #
    # only run if the output VCF doesn't already exist
    if not AnyPath(REHEADERED_OUT).exists():
        logging.info(
            f'The Reheadered VCF "{REHEADERED_OUT}" doesn\'t exist; regenerating'
        )

        if not AnyPath(HAIL_VCF_OUT).exists():
            logging.info(
                f'The Labelled VCF "{HAIL_VCF_OUT}" doesn\'t exist; regenerating'
            )

            # do we need to run the full annotation stage?
            if not AnyPath(matrix_path.rstrip('/') + '/').exists():
                raise Exception(
                    f'Currently this process demands an annotated '
                    f'MatrixTable. The provided path "{matrix_path}" '
                    f'does not exist or is inaccessible'
                )

            prior_job = handle_hail_filtering(
                batch=batch,
                matrix_path=matrix_path,
                config=config_json,
                prior_job=prior_job,
            )

        # --------------------------------- #
        # bcftools re-headering of hail VCF #
        # --------------------------------- #
        # copy the labelled output file into the remaining batch jobs
        hail_output_in_batch = batch.read_input_group(
            **{'vcf.bgz': HAIL_VCF_OUT, 'vcf.bgz.tbi': HAIL_VCF_OUT + '.tbi'}
        )

        # this is no longer explicitly required...
        # it was required to run slivar: geneId, consequences, and transcript
        # keeping in to ensure the VCF can be interpreted without the config
        bcftools_job = handle_reheader_job(
            batch=batch,
            local_vcf=hail_output_in_batch['vcf.bgz'],
            config_dict=config_dict,
            prior_job=prior_job,
        )
        batch.write_output(bcftools_job.vcf, REHEADERED_PREFIX)
        reheadered_vcf_in_batch = bcftools_job.vcf['vcf.bgz']

    # if it exists remotely, read into a batch
    else:
        vcf_in_batch = batch.read_input_group(
            **{'vcf.bgz': HAIL_VCF_OUT, 'vcf.bgz.tbi': HAIL_VCF_OUT + '.tbi'}
        )
        reheadered_vcf_in_batch = vcf_in_batch['vcf.bgz']

    # use compound-hets and labelled VCF to identify plausibly pathogenic
    # variants where the MOI is viable compared to the PanelApp expectation
    _results_job = handle_results_job(
        batch=batch,
        config=config_json,
        comp_het=COMP_HET_JSON,
        reheadered_vcf=reheadered_vcf_in_batch,
        pedigree=pedigree,
        prior_job=prior_job,
    )
    batch.run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )
    main()  # pylint: disable=E1120
