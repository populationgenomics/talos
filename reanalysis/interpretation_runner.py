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
"""
from datetime import datetime
import json
import logging
import os
import sys

import click
from cloudpathlib import AnyPath, CloudPath
import hailtop.batch as hb

from cpg_utils.config import get_config
from cpg_utils.deploy_config import get_deploy_config
from cpg_utils.git import (
    prepare_git_job,
    get_git_commit_ref_of_current_repository,
    get_organisation_name_from_current_directory,
    get_repo_name_from_current_directory,
)
from cpg_utils.hail_batch import (
    authenticate_cloud_credentials_in_job,
    copy_common_env,
    output_path,
    query_command,
    remote_tmpdir,
    image_path,
)

import annotation
from utils import read_json_from_path, FileTypes, identify_file_type
from vep.jobs import vep_jobs, SequencingType

# BIG TODO, HACK HERE
def fix_output_path(input: str) -> str:
    return input.replace('/severalgenomes', '/cpg-severalgenomes', 1)
    
# static paths to write outputs
INPUT_AS_VCF = fix_output_path(output_path('prior_to_annotation.vcf.bgz'))

# phases of annotation
ANNOTATED_MT = fix_output_path(output_path('annotated_variants.mt'))

# panelapp query results
PANELAPP_JSON_OUT = fix_output_path(output_path('panelapp_data.json'))

# output of labelling task in Hail
HAIL_VCF_OUT = fix_output_path(output_path('hail_categorised.vcf.bgz'))


# local script references
HAIL_FILTER = os.path.join(os.path.dirname(__file__), 'hail_filter_and_label.py')
HTML_SCRIPT = os.path.join(os.path.dirname(__file__), 'html_builder.py')
QUERY_PANELAPP = os.path.join(os.path.dirname(__file__), 'query_panelapp.py')
RESULTS_SCRIPT = os.path.join(os.path.dirname(__file__), 'validate_categories.py')
MT_TO_VCF_SCRIPT = os.path.join(os.path.dirname(__file__), 'mt_to_vcf.py')


def set_job_resources(
    job: hb.batch.job.Job,
    auth=False,
    git=False,
    prior_job: hb.batch.job.Job | None = None,
    memory: str = 'standard',
):
    """
    applied resources to the job
    :param job:
    :param auth: if true, authenticate gcloud in this container
    :param git: if true, pull this repository into container
    :param prior_job:
    :param memory:
    """
    # apply all settings
    job.cpu(2).image(image_path('hail')).memory(memory).storage('20G')

    if prior_job is not None:
        job.depends_on(prior_job)

    if auth:
        # TODO, the cloud check should really be done in CPG utils, also there's a question of whether this is actually necessary - I think it's done
        # in GCP to handle use of AnyPath. We may not actually need this because of how HailAzureCloudPath handles auth.
        cloud = get_deploy_config().to_dict()['cloud']
        if cloud == 'gcp':
            authenticate_cloud_credentials_in_job(job)
        elif cloud == 'azure':
            raise NotImplementedError()

    if git:
        # copy the relevant scripts into a Driver container instance
        prepare_git_job(
            job=job,
            organisation=get_organisation_name_from_current_directory(),
            repo_name=get_repo_name_from_current_directory(),
            commit=get_git_commit_ref_of_current_repository(),
        )


def mt_to_vcf(batch: hb.Batch, input_file: str):
    """
    takes a MT and converts to VCF
    :param batch:
    :param input_file:
    :return:
    """
    mt_to_vcf_job = batch.new_job(name='Convert MT to VCF')
    auth = get_deploy_config().to_dict()['cloud'] == 'gcp'
    set_job_resources(mt_to_vcf_job, git=True, auth=auth)

    job_cmd = (
        f'PYTHONPATH=$(pwd) python3 {MT_TO_VCF_SCRIPT} '
        f'--input {input_file} '
        f'--output {INPUT_AS_VCF}'
    )

    logging.info(f'Command used to convert MT: {job_cmd}')
    copy_common_env(mt_to_vcf_job)
    mt_to_vcf_job.command(job_cmd)
    return mt_to_vcf_job


def annotate_vcf(
    input_vcf: str,
    batch: hb.Batch,
    vep_temp: str,
    vep_out: str,
    seq_type: SequencingType = SequencingType.GENOME,
) -> list[hb.batch.job.Job]:
    """
    takes the VCF path, schedules all annotation jobs, creates MT with VEP annos.

    should this be separated out into a script and run end2end, or should we
    continue in this same runtime? These jobs are scheduled into this batch with
    appropriate dependencies, so keeping this structure seems valid

    :param input_vcf:
    :param batch:
    :param vep_temp:
    :param vep_out:
    :param seq_type:
    :return:
    """

    # generate the jobs which run VEP & collect the results
    return vep_jobs(
        b=batch,
        vcf_path=AnyPath(input_vcf),
        tmp_bucket=AnyPath(vep_temp),
        out_path=AnyPath(vep_out),
        overwrite=False,  # don't re-run annotation on completed chunks
        sequencing_type=seq_type,
        job_attrs={},
    )


def annotated_mt_from_ht_and_vcf(
    input_vcf: str,
    batch: hb.Batch,
    vep_ht: str,
    job_attrs: dict | None = None,
) -> hb.batch.job.Job:
    """
    apply the HT of annotations to the VCF, save as MT
    :return:
    """
    apply_anno_job = batch.new_job('HT + VCF = MT', job_attrs)

    copy_common_env(apply_anno_job)
    apply_anno_job.image(image_path('hail'))

    cmd = query_command(
        annotation,
        annotation.apply_annotations.__name__,
        input_vcf,
        vep_ht,
        ANNOTATED_MT,
        setup_gcp=True,
        packages=['seqr-loader==1.2.5'],
    )
    apply_anno_job.command(cmd)
    return apply_anno_job


def handle_panelapp_job(
    batch: hb.Batch,
    extra_panel: tuple[str],
    gene_list: str | None = None,
    prior_job: hb.batch.job.Job | None = None,
) -> hb.batch.job.Job:
    """

    :param batch:
    :param extra_panel:
    :param gene_list:
    :param prior_job:
    """
    panelapp_job = batch.new_job(name='query panelapp')
    auth = get_deploy_config().to_dict()['cloud'] == 'gcp'
    set_job_resources(panelapp_job, auth=auth, git=True, prior_job=prior_job)

    panelapp_command = f'python3 {QUERY_PANELAPP} --out_path {PANELAPP_JSON_OUT} '
    if gene_list is not None:
        panelapp_command += f'--gene_list {gene_list} '
    if extra_panel is not None and len(extra_panel) != 0:
        panelapp_command += f'-p {" ".join(extra_panel)} '

    if prior_job is not None:
        panelapp_job.depends_on(prior_job)

    logging.info(f'PanelApp Command: {panelapp_command}')
    panelapp_job.command(panelapp_command)
    return panelapp_job


def handle_hail_filtering(
    batch: hb.Batch,
    config: str,
    plink_file: str,
    prior_job: hb.batch.job.Job | None = None,
) -> hb.batch.job.BashJob:
    """
    hail-query backend version of the filtering implementation
    use the init query service instead of running inside dataproc

    :param batch:
    :param config:
    :param plink_file:
    :param prior_job:
    :return:
    """

    labelling_job = batch.new_job(name='hail filtering')
    auth = get_deploy_config().to_dict()['cloud'] == 'gcp'
    set_job_resources(
        labelling_job, auth=auth, git=True, prior_job=prior_job, memory='16Gi'
    )
    labelling_command = (
        f'pip install . && '
        f'python3 {HAIL_FILTER} '
        f'--mt {ANNOTATED_MT} '
        f'--panelapp {PANELAPP_JSON_OUT} '
        f'--config_path {config} '
        f'--plink {plink_file}'
    )

    logging.info(f'Labelling Command: {labelling_command}')
    labelling_job.command(labelling_command)
    copy_common_env(labelling_job)
    return labelling_job


def handle_results_job(
    batch: hb.Batch,
    config: str,
    labelled_vcf: str,
    pedigree: str,
    output_dict: dict[str, dict[str, str]],
    prior_job: hb.batch.job.Job | None = None,
) -> hb.batch.job.Job:
    """
    one container to run the MOI checks, and the presentation

    :param batch:
    :param config:
    :param labelled_vcf:
    :param pedigree:
    :param output_dict: paths to the
    :param prior_job:
    :return:
    """

    results_job = batch.new_job(name='finalise_results')
    auth = get_deploy_config().to_dict()['cloud'] == 'gcp'
    set_job_resources(results_job, auth=True, git=True, prior_job=prior_job)
    results_command = (
        'pip install . && '
        f'python3 {RESULTS_SCRIPT} '
        f'--config_path {config} '
        f'--labelled_vcf {labelled_vcf} '
        f'--panelapp {PANELAPP_JSON_OUT} '
        f'--pedigree {pedigree} '
        f'--out_json {output_dict["results"]} && '
        f'python3 {HTML_SCRIPT} '
        f'--results {output_dict["results"]} '
        f'--config_path {config} '
        f'--panelapp {PANELAPP_JSON_OUT} '
        f'--pedigree {pedigree} '
        f'--out_path {output_dict["web_html"]}'
    )
    logging.info(f'Results command: {results_command}')
    results_job.command(results_command)
    return results_job


# def cpg_utils_resolve_blob_name_as_hail_path(blob_name: str) -> str:
#     account = 'sevgen002'
#     container = 'cpg-severalgenomes-main'
#     return f'hail-az://{account}/{container}/{blob_name.lstrip("/")}'


@click.command()
@click.option(
    '--input_path', help='variant matrix table or VCF to analyse', required=True
)
@click.option('--config_json', help='JSON dict of runtime settings', required=True)
@click.option('--plink_file', help='Plink file path for the cohort', required=True)
@click.option(
    '--extra_panel',
    help='Any additional panelapp IDs to add to the Mendeliome. '
    'Multiple can be added as "--extra_panel 123 --extra_panel 456',
    required=False,
    multiple=True,
)
@click.option(
    '--panel_genes', help='JSON Gene list for use in analysis', required=False
)
@click.option(
    '--singletons', help='location of a plink file for the singletons', required=False
)
@click.option(
    '--skip_annotation',
    help='if set, a MT with appropriate annotations can be provided',
    is_flag=True,
    default=False,
)
def main(
    input_path: str,
    config_json: str,
    plink_file: str,
    extra_panel: tuple[str],
    panel_genes: str | None = None,
    singletons: str | None = None,
    skip_annotation: bool = False,
):
    """
    main method, which runs the full reanalysis process

    :param input_path: annotated input matrix table or VCF
    :param config_json:
    :param plink_file:
    :param panel_genes:
    :param extra_panel:
    :param singletons:
    :param skip_annotation:
    """

    if not AnyPath(input_path).exists():
        raise Exception(
            f'The provided path "{input_path}" does not exist or is inaccessible'
        )

    logging.info('Starting the reanalysis batch')

    config_dict = read_json_from_path(config_json)
    config_dict.update(
        {
            'latest_run': f'{datetime.now():%Y-%m-%d %H:%M:%S%z}',
            'input_file': input_path,
            'panelapp_file': PANELAPP_JSON_OUT,
            'cohort': get_config()['workflow']['dataset'],
        }
    )

    # create output paths with optional suffixes
    vep_stage_tmp = output_path('vep_temp', config_dict.get('tmp_suffix') or None)
    vep_ht_tmp = output_path(
        'vep_annotations.ht', config_dict.get('tmp_suffix') or None
    )

    # separate paths for familial and singleton analysis
    output_dict = {
        'default': {
            'web_html': output_path(
                'summary_output.html', config_dict.get('web_suffix') or None
            ),
            'results': output_path('summary_results.json'),
        },
        'singletons': {
            'web_html': output_path(
                'singleton_output.html', config_dict.get('web_suffix') or None
            ),
            'results': output_path('singleton_results.json'),
        },
    }

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

    # read the ped file into the Batch
    pedigree_in_batch = batch.read_input(plink_file)

    # set a first job in this batch
    prior_job = None

    # -------------------------- #
    # Convert MT to a VCF format #
    # -------------------------- #
    # determine the input type - if MT, decompose to VCF prior to annotation
    input_file_type = identify_file_type(input_path)
    assert input_file_type in [
        FileTypes.VCF_GZ,
        FileTypes.VCF_BGZ,
        FileTypes.MATRIX_TABLE,
    ], (
        f'inappropriate input type provided: {input_file_type}; '
        f'this is designed for MT or compressed VCF only'
    )

    if input_file_type == FileTypes.MATRIX_TABLE:
        if skip_annotation:
            config_dict.update({'aip_annotated': False})
            # overwrite the expected annotation output path
            global ANNOTATED_MT  # pylint: disable=W0603
            ANNOTATED_MT = input_path

        else:
            prior_job = mt_to_vcf(batch=batch, input_file=input_path)
            config_dict.update({'vcf_created': INPUT_AS_VCF})
            # overwrite input path with file we just created
            input_path = INPUT_AS_VCF

    # ------------------------------------- #
    # split the VCF, and annotate using VEP #
    # ------------------------------------- #
    if not CloudPath(ANNOTATED_MT).exists():
        # need to run the annotation phase
        # uses default values from RefData
        annotation_jobs = annotate_vcf(
            input_path, batch=batch, vep_temp=vep_stage_tmp, vep_out=vep_ht_tmp
        )

        # if convert-to-VCF job exists, assign as an annotation dependency
        if prior_job:
            for job in annotation_jobs:
                job.depends_on(prior_job)

        # apply annotations
        prior_job = annotated_mt_from_ht_and_vcf(
            input_vcf=input_path, batch=batch, job_attrs={}, vep_ht=vep_ht_tmp
        )
        prior_job.depends_on(*annotation_jobs)

        config_dict.update({'aip_annotated': True})
    else:
        logging.info("Using previously-annotated MT")

    # -------------------------------- #
    # query panelapp for panel details #
    # -------------------------------- #

    if not AnyPath(PANELAPP_JSON_OUT).exists():
        logging.info(f"PanelApp JSON {PANELAPP_JSON_OUT} doesn\'t exist, generating.")
        prior_job = handle_panelapp_job(
            batch=batch,
            extra_panel=extra_panel,
            gene_list=panel_genes,
            prior_job=prior_job,
        )
    else:
        logging.info("Using previous PanelApp JSON")

    # ----------------------- #
    # run hail categorisation #
    # ----------------------- #
    if not AnyPath(HAIL_VCF_OUT).exists():
        logging.info(f'The Labelled VCF "{HAIL_VCF_OUT}" doesn\'t exist; regenerating')
        prior_job = handle_hail_filtering(
            batch=batch,
            config=config_json,
            prior_job=prior_job,
            plink_file=pedigree_in_batch,
        )
    else:
        logging.info(f"Using previous labelled VCF: {HAIL_VCF_OUT}")

    # read that VCF into the batch as a local file
    labelled_vcf_in_batch = batch.read_input_group(
        vcf=HAIL_VCF_OUT, tbi=HAIL_VCF_OUT + '.tbi'
    ).vcf

    # if singleton PED supplied, also run as singletons w/separate outputs
    analysis_rounds = [(pedigree_in_batch, 'default')]
    if singletons and AnyPath(singletons).exists():
        pedigree_singletons = batch.read_input(singletons)
        analysis_rounds.append((pedigree_singletons, 'singletons'))
    else:
        logging.info("Skipping singleton analysis")

    # pointing this analysis at the updated config file, including input metadata
    for relationships, analysis_index in analysis_rounds:
        logging.info(f'running analysis in {analysis_index} mode')
        _results_job = handle_results_job(
            batch=batch,
            config=output_path('latest_config.json'),
            labelled_vcf=labelled_vcf_in_batch,
            pedigree=relationships,
            output_dict=output_dict[analysis_index],
            prior_job=prior_job,
        )

    # save the json file into the batch output, with latest run details
    with AnyPath(output_path('latest_config.json')).open('w') as handle:
        json.dump(config_dict, handle , indent=True)

    # write pedigree content to the output folder
    with AnyPath(output_path('latest_pedigree.fam')).open('w') as handle:
        handle.writelines(AnyPath(plink_file).open().readlines())

    if singletons:
        with AnyPath(output_path('latest_singletons.fam')).open('w') as handle:
            handle.writelines(AnyPath(singletons).open().readlines())

    batch.run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )
    main()  # pylint: disable=E1120
