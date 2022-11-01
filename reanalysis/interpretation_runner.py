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


from argparse import ArgumentParser
from datetime import datetime
import logging
import os
import sys

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config
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
from utils import FileTypes, identify_file_type
from vep.jobs import vep_jobs, SequencingType

# exact time that this run occurred
EXECUTION_TIME = f'{datetime.now():%Y-%m-%d %H:%M}'

# static paths to write outputs
INPUT_AS_VCF = output_path('prior_to_annotation.vcf.bgz')

# phases of annotation
ANNOTATED_MT = output_path('annotated_variants.mt')

# panelapp query results
PANELAPP_JSON_OUT = output_path(
    'panelapp_data', get_config()['buckets'].get('analysis_suffix')
)

# output of labelling task in Hail
HAIL_VCF_OUT = output_path('hail_categorised.vcf.bgz')


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

    # copy the env variables into the container
    # specifically the CPG_CONFIG_PATH value
    copy_common_env(job)

    if prior_job is not None:
        job.depends_on(prior_job)

    if auth:
        authenticate_cloud_credentials_in_job(job)

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
    set_job_resources(mt_to_vcf_job, git=True, auth=True)

    job_cmd = (
        f'PYTHONPATH=$(pwd) python3 {MT_TO_VCF_SCRIPT} '
        f'--input {input_file} '
        f'--output {INPUT_AS_VCF}'
    )

    logging.info(f'Command used to convert MT: {job_cmd}')
    mt_to_vcf_job.command(job_cmd)
    return mt_to_vcf_job


def annotate_vcf(
    input_vcf: str,
    batch: hb.Batch,
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
    :param vep_out:
    :param seq_type:
    :return:
    """

    # generate the jobs which run VEP & collect the results
    return vep_jobs(
        b=batch,
        vcf_path=to_path(input_vcf),
        tmp_bucket=to_path(
            output_path('vep_temp', get_config()['buckets'].get('tmp_suffix'))
        ),
        out_path=to_path(vep_out),
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
    extra_panels: list[str] | None = None,
    participant_panels: str | None = None,
    prior_job: hb.batch.job.Job | None = None,
) -> hb.batch.job.Job:
    """

    :param batch:
    :param extra_panels:
    :param participant_panels:
    :param prior_job:
    """
    panelapp_job = batch.new_job(name='query panelapp')
    set_job_resources(panelapp_job, auth=True, git=True, prior_job=prior_job)

    panelapp_command = f'python3 {QUERY_PANELAPP} --out_path {PANELAPP_JSON_OUT} '

    if extra_panels:
        panelapp_command += f'-p {" ".join(extra_panels)} '

    if participant_panels:
        panelapp_command += f'--panel_file {participant_panels} '

    if prior_job is not None:
        panelapp_job.depends_on(prior_job)

    logging.info(f'PanelApp Command: {panelapp_command}')
    panelapp_job.command(panelapp_command)
    return panelapp_job


def handle_hail_filtering(
    batch: hb.Batch,
    plink_file: str,
    prior_job: hb.batch.job.Job | None = None,
) -> hb.batch.job.BashJob:
    """
    hail-query backend version of the filtering implementation
    use the init query service instead of running inside dataproc

    :param batch:
    :param plink_file:
    :param prior_job:
    :return:
    """

    labelling_job = batch.new_job(name='hail filtering')
    set_job_resources(
        labelling_job, auth=True, git=True, prior_job=prior_job, memory='16Gi'
    )
    labelling_command = (
        f'pip install . && '
        f'python3 {HAIL_FILTER} '
        f'--mt {ANNOTATED_MT} '
        f'--panelapp {PANELAPP_JSON_OUT}.json '
        f'--plink {plink_file}'
    )

    logging.info(f'Labelling Command: {labelling_command}')
    labelling_job.command(labelling_command)
    return labelling_job


def handle_results_job(
    batch: hb.Batch,
    labelled_vcf: str,
    pedigree: str,
    output_dict: dict[str, dict[str, str]],
    prior_job: hb.batch.job.Job | None = None,
    participant_panels: str | None = None,
    input_path: str | None = None,
) -> hb.batch.job.Job:
    """
    one container to run the MOI checks, and the presentation

    :param batch:
    :param labelled_vcf:
    :param pedigree:
    :param output_dict: paths to the
    :param prior_job:
    :param participant_panels: JSON of relevant panels per participant
    :param input_path: source file for the analysis process
    :return:
    """

    results_job = batch.new_job(name='finalise_results')
    set_job_resources(results_job, auth=True, git=True, prior_job=prior_job)

    gene_filter_files = (
        (
            f'--participant_panels {participant_panels} '
            f'--panel_genes {PANELAPP_JSON_OUT}_per_panel.json '
        )
        if participant_panels
        else ''
    )

    report_from_file = (
        f'--results {output_dict["results"]}_panel_filtered.json '
        if participant_panels
        else f'--results {output_dict["results"]}_full.json '
    )

    path_input = f'--input_path {input_path} ' if input_path else ''

    results_command = (
        'pip install . && '
        f'python3 {RESULTS_SCRIPT} '
        f'--labelled_vcf {labelled_vcf} '
        f'--panelapp {PANELAPP_JSON_OUT}.json '
        f'--pedigree {pedigree} '
        f'--out_json {output_dict["results"]} '
        f'{gene_filter_files} {path_input} && '
        f'python3 {HTML_SCRIPT} '
        f'{report_from_file} '
        f'--panelapp {PANELAPP_JSON_OUT}.json '
        f'--pedigree {pedigree} '
        f'--out_path {output_dict["web_html"]}'
    )
    logging.info(f'Results command: {results_command}')
    results_job.command(results_command)
    return results_job


def main(
    input_path: str,
    pedigree: str,
    extra_panels: list[str],
    participant_panels: str | None,
    singletons: str | None = None,
    skip_annotation: bool = False,
):
    """
    main method, which runs the full reanalysis process

    :param input_path: annotated input matrix table or VCF
    :param pedigree:
    :param extra_panels:
    :param participant_panels:
    :param singletons:
    :param skip_annotation:
    """

    assert to_path(
        input_path
    ).exists(), f'The provided path "{input_path}" does not exist or is inaccessible'

    logging.info('Starting the reanalysis batch')

    # region: output files lookup
    # separate paths for familial and singleton analysis
    output_dict = {
        'default': {
            'web_html': output_path(
                'summary_output.html', get_config()['buckets'].get('web_suffix')
            ),
            'results': output_path(
                'summary_results', get_config()['buckets'].get('analysis_suffix')
            ),
        },
        'singletons': {
            'web_html': output_path(
                'singleton_output.html', get_config()['buckets'].get('web_suffix')
            ),
            'results': output_path(
                'singleton_results', get_config()['buckets'].get('analysis_suffix')
            ),
        },
    }
    # endregion

    # region: batch instantiation
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
    # endregion

    # read the ped file into the Batch
    pedigree_in_batch = batch.read_input(pedigree)

    # set a first job in this batch
    prior_job = None

    # region: MT to VCF
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
            # overwrite the expected annotation output path
            global ANNOTATED_MT  # pylint: disable=W0603
            ANNOTATED_MT = input_path

        else:
            prior_job = mt_to_vcf(batch=batch, input_file=input_path)
            # overwrite input path with file we just created
            input_path = INPUT_AS_VCF
    # endregion

    # region: split & annotate VCF
    if not to_path(ANNOTATED_MT).exists():
        # need to run the annotation phase
        # uses default values from RefData

        vep_ht_tmp = output_path(
            'vep_annotations.ht', get_config()['buckets'].get('tmp_suffix')
        )
        annotation_jobs = annotate_vcf(input_path, batch=batch, vep_out=vep_ht_tmp)

        # if convert-to-VCF job exists, assign as an annotation dependency
        if prior_job:
            for job in annotation_jobs:
                job.depends_on(prior_job)

        # apply annotations
        prior_job = annotated_mt_from_ht_and_vcf(
            input_vcf=input_path, batch=batch, job_attrs={}, vep_ht=vep_ht_tmp
        )
        prior_job.depends_on(*annotation_jobs)

    # endregion

    #  region: query panelapp
    if (not to_path(f'{PANELAPP_JSON_OUT}.json').exists()) or (
        participant_panels
        and not to_path(f'{PANELAPP_JSON_OUT}_per_panel.json').exists()
    ):
        prior_job = handle_panelapp_job(
            batch=batch,
            extra_panels=extra_panels,
            participant_panels=participant_panels,
            prior_job=prior_job,
        )
    # endregion

    # region: hail categorisation
    if not to_path(HAIL_VCF_OUT).exists():
        logging.info(f'The Labelled VCF "{HAIL_VCF_OUT}" doesn\'t exist; regenerating')
        prior_job = handle_hail_filtering(
            batch=batch,
            prior_job=prior_job,
            plink_file=pedigree_in_batch,
        )
    # endregion

    # read VCF into the batch as a local file
    labelled_vcf_in_batch = batch.read_input_group(
        vcf=HAIL_VCF_OUT, tbi=HAIL_VCF_OUT + '.tbi'
    ).vcf

    # region: singleton decisions
    # if singleton PED supplied, also run as singletons w/separate outputs
    analysis_rounds = [(pedigree_in_batch, 'default')]
    if singletons and to_path(singletons).exists():
        to_path(singletons).copy(
            output_path(
                f'singletons_{EXECUTION_TIME}.fam',
                get_config()['buckets'].get('analysis_suffix'),
            )
        )
        pedigree_singletons = batch.read_input(singletons)
        analysis_rounds.append((pedigree_singletons, 'singletons'))
    # endregion

    # region: run results job
    # pointing this analysis at the updated config file, including input metadata
    for relationships, analysis_index in analysis_rounds:
        logging.info(f'running analysis in {analysis_index} mode')
        _results_job = handle_results_job(
            batch=batch,
            labelled_vcf=labelled_vcf_in_batch,
            pedigree=relationships,
            output_dict=output_dict[analysis_index],
            prior_job=prior_job,
            participant_panels=participant_panels,
            input_path=input_path,
        )
    # endregion

    # region: copy data out
    # if we ran with per-participant panel data, copy to output folder
    # include datetime to differentiate output files and prevent clashes
    if participant_panels:
        to_path(participant_panels).copy(
            output_path(
                f'pid_to_panels_{EXECUTION_TIME}.json',
                get_config()['buckets'].get('analysis_suffix'),
            )
        )

    # write pedigree content to the output folder
    to_path(pedigree).copy(
        output_path(
            f'pedigree_{EXECUTION_TIME}.fam',
            get_config()['buckets'].get('analysis_suffix'),
        )
    )
    # endregion

    batch.run(wait=False)


if __name__ == '__main__':
    import debugpy

    debugpy.listen(('0.0.0.0', 5678))
    debugpy.wait_for_client()

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )

    parser = ArgumentParser()
    parser.add_argument('-i', help='variant data to analyse', required=True)
    parser.add_argument('--pedigree', help='in Plink format', required=True)
    parser.add_argument('--singletons', help='singletons in Plink format')
    panel_args = parser.add_mutually_exclusive_group()
    panel_args.add_argument(
        '--extra_panels', help='any additional panel IDs', nargs='+', default=[]
    )
    panel_args.add_argument(
        '--participant_panels',
        help='JSON file containing per-participant panel details',
    )
    parser.add_argument(
        '--skip_annotation',
        help='if set, annotation will not be repeated',
        action='store_true',
    )
    args = parser.parse_args()
    main(
        input_path=args.i,
        pedigree=args.pedigree,
        extra_panels=args.extra_panels,
        participant_panels=args.participant_panels,
        skip_annotation=args.skip_annotation,
    )
