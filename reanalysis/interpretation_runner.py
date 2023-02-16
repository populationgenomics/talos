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


# pylint: disable=too-many-branches


import logging
import os
import sys
from argparse import ArgumentParser
from datetime import datetime

from hailtop.batch.job import BashJob, Job

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import (
    authenticate_cloud_credentials_in_job,
    copy_common_env,
    dataset_path,
    output_path,
)
from cpg_utils.git import get_git_root_relative_path_from_absolute

from cpg_workflows.batch import get_batch
from cpg_workflows.jobs.seqr_loader import annotate_cohort_jobs
from cpg_workflows.jobs.vep import add_vep_jobs
from cpg_workflows.jobs.joint_genotyping import add_make_sitesonly_job

from reanalysis import (
    hail_filter_and_label,
    html_builder,
    mt_to_vcf,
    query_panelapp,
    summarise_clinvar_entries,
    validate_categories,
)
from reanalysis.utils import FileTypes, identify_file_type


# exact time that this run occurred
EXECUTION_TIME = f'{datetime.now():%Y-%m-%d_%H:%M}'

# static paths to write outputs
INPUT_AS_VCF = output_path('prior_to_annotation.vcf.bgz')

# phases of annotation
ANNOTATED_MT = output_path('annotated_variants.mt')

# panelapp query results
PANELAPP_JSON_OUT = output_path('panelapp_data.json', 'analysis')

# output of labelling task in Hail
HAIL_VCF_OUT = output_path('hail_categorised.vcf.bgz', 'analysis')


def set_job_resources(job: Job, prior_job: Job | None = None, memory: str = 'standard'):
    """
    applied resources to the job

    Args:
        job ():
        prior_job ():
        memory (str): lowmem/standard/highmem
    """
    # apply all settings
    job.cpu(2).image(get_config()['workflow']['driver_image']).memory(memory).storage(
        '20G'
    )

    # copy the env variables into the container; specifically CPG_CONFIG_PATH
    copy_common_env(job)

    if prior_job is not None:
        if isinstance(prior_job, list):
            job.depends_on(*prior_job)
        else:
            job.depends_on(prior_job)

    authenticate_cloud_credentials_in_job(job)


def handle_clinvar() -> tuple[Job | None, str]:
    """
    set up a job to handle the clinvar summarising

    Returns:
        the batch job for creating the new summary
        the path to the current clinvar summary file
    """

    # is it time to re-process clinvar?
    clinvar_prefix = dataset_path(
        f'clinvar_summaries/{datetime.now().strftime("%Y_%m")}'
    )
    clinvar_summary = os.path.join(clinvar_prefix, 'clinvar.ht')

    if to_path(clinvar_summary).exists():
        return None, clinvar_summary

    # create a bash job to copy data from remote
    bash_job = get_batch().new_bash_job(name='copy clinvar files to local')
    set_job_resources(bash_job)

    directory = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/'
    sub_file = 'submission_summary.txt.gz'
    var_file = 'variant_summary.txt.gz'

    bash_job.command(
        (
            f'wget -q {directory}{sub_file} -O {bash_job.subs} && '
            f'wget -q {directory}{var_file} -O {bash_job.vars}'
        )
    )

    # write output files
    get_batch().write_output(bash_job.subs, os.path.join(clinvar_prefix, sub_file))
    get_batch().write_output(bash_job.vars, os.path.join(clinvar_prefix, var_file))

    # create a job to run the summary
    summarise = get_batch().new_job(name='summarise clinvar')
    set_job_resources(summarise, prior_job=bash_job)

    script_path = get_git_root_relative_path_from_absolute(summarise_clinvar_entries.__file__)
    summarise.command(
        f'python3 {script_path} '
        f'-s {bash_job.subs} '
        f'-v {bash_job.vars} '
        f'-o {os.path.join(clinvar_prefix, "clinvar.ht")}'
    )
    return summarise, clinvar_summary


def setup_mt_to_vcf(input_file: str) -> Job:
    """
    set up a job MatrixTable conversion to VCF, prior to annotation

    Args:
        input_file (str): path to the MatrixTable

    Returns:
        the new job, available for dependency setting
    """

    job = get_batch().new_job(name='Convert MT to VCF')
    set_job_resources(job)

    script_path = get_git_root_relative_path_from_absolute(mt_to_vcf.__file__)
    cmd = f'python3 {script_path} --input {input_file} --output {INPUT_AS_VCF}'

    logging.info(f'Command used to convert MT: {cmd}')
    job.command(cmd)
    return job


def handle_panelapp_job(
    participant_panels: str | None = None,
    previous: str | None = None,
    prior_job: Job | None = None,
) -> Job:
    """
    creates and runs the panelapp query job

    Args:
        participant_panels (str):
        previous (str): optional, path to a prior gene list
        prior_job ():

    Returns:
        the Job, which other parts of the workflow may become dependent on
    """
    panelapp_job = get_batch().new_job(name='query panelapp')
    set_job_resources(panelapp_job, prior_job=prior_job)
    copy_common_env(panelapp_job)

    script_path = get_git_root_relative_path_from_absolute(query_panelapp.__file__)
    query_cmd = f'python3 {script_path} --out_path {PANELAPP_JSON_OUT} '

    if participant_panels is not None:
        query_cmd += f'--panels {participant_panels} '

    if previous is not None:
        query_cmd += f'--previous {previous} '

    logging.info(f'PanelApp Command: {query_cmd}')
    panelapp_job.command(query_cmd)
    return panelapp_job


def handle_hail_filtering(
    plink_file: str, clinvar: str, prior_job: Job | None = None
) -> BashJob:
    """
    hail-query backend version of the filtering implementation
    use the init query service instead of running inside dataproc

    Args:
        plink_file (str): path to a pedigree
        clinvar (str): path to the clinvar re-summary
        prior_job ():

    Returns:
        the Batch job running the hail filtering process
    """

    labelling_job = get_batch().new_job(name='hail filtering')
    set_job_resources(labelling_job, prior_job=prior_job, memory='16Gi')
    script_path = get_git_root_relative_path_from_absolute(hail_filter_and_label.__file__)
    labelling_command = (
        f'python3 {script_path} '
        f'--mt {ANNOTATED_MT} '
        f'--panelapp {PANELAPP_JSON_OUT} '
        f'--plink {plink_file} '
        f'--clinvar {clinvar} '
    )

    logging.info(f'Labelling Command: {labelling_command}')
    labelling_job.command(labelling_command)
    return labelling_job


def handle_results_job(
    labelled_vcf: str,
    pedigree: str,
    input_path: str,
    output_dict: dict[str, dict[str, str]],
    prior_job: Job | None = None,
    participant_panels: str | None = None,
):
    """
    one container to run the MOI checks, and the presentation

    Args:
        labelled_vcf ():
        pedigree ():
        input_path (str): path to the input file, logged in metadata
        output_dict ():
        prior_job ():
        participant_panels ():
    """

    results_job = get_batch().new_job(name='finalise_results')
    set_job_resources(results_job, prior_job=prior_job)

    gene_filter_files = (
        f'--participant_panels {participant_panels} ' if participant_panels else ''
    )

    validation_script_path = get_git_root_relative_path_from_absolute(validate_categories.__file__)
    html_script_path = get_git_root_relative_path_from_absolute(html_builder.__file__)

    results_command = (
        f'python3 {validation_script_path} '
        f'--labelled_vcf {labelled_vcf} '
        f'--panelapp {PANELAPP_JSON_OUT} '
        f'--pedigree {pedigree} '
        f'--out_json {output_dict["results"]} '
        f'--input_path {input_path} '
        f'{gene_filter_files} && '
        f'python3 {html_script_path} '
        f'--results {output_dict["results"]} '
        f'--panelapp {PANELAPP_JSON_OUT} '
        f'--pedigree {pedigree} '
        f'--out_path {output_dict["web_html"]}'
    )
    logging.info(f'Results command: {results_command}')
    results_job.command(results_command)


def main(
    input_path: str,
    pedigree: str,
    participant_panels: str | None,
    previous: str | None,
    singletons: str | None = None,
    skip_annotation: bool = False,
):
    """
    main method, which runs the full reanalysis process

    Args:
        input_path (): path to the VCF/MT
        pedigree (): family file for this analysis
        participant_panels (): file containing panels-per-family (optional)
        previous (): gene panel data from prior analysis (optional)
        singletons (): optional second Pedigree file without families
        skip_annotation (): if the input is annotated, don't re-run
    """

    assert to_path(
        input_path
    ).exists(), f'The provided path {input_path!r} does not exist or is inaccessible'

    logging.info('Starting the reanalysis batch')

    # region: output files lookup
    # separate paths for familial and singleton analysis
    output_dict = {
        'default': {
            'web_html': output_path('summary_output.html', 'web'),
            'results': output_path('summary_results.json', 'analysis'),
        },
        'singletons': {
            'web_html': output_path('singleton_output.html', 'web'),
            'results': output_path('singleton_results.json', 'analysis'),
        },
    }
    # endregion

    # find clinvar table, and re-process if required
    prior_job, clinvar_table = handle_clinvar()

    # region: MT to VCF
    # determine the input type - if MT, decompose to VCF prior to annotation
    input_file_type = identify_file_type(input_path)
    assert input_file_type in [
        FileTypes.VCF_GZ,
        FileTypes.VCF_BGZ,
        FileTypes.MATRIX_TABLE,
    ], (
        f'inappropriate input type provided: {input_file_type}; '
        'this is designed for MT or compressed VCF only'
    )

    if input_file_type == FileTypes.MATRIX_TABLE:
        if skip_annotation:
            # overwrite the expected annotation output path
            global ANNOTATED_MT  # pylint: disable=W0603
            ANNOTATED_MT = input_path

        else:
            if not to_path(INPUT_AS_VCF).exists():
                prior_job = setup_mt_to_vcf(input_file=input_path)

            # overwrite input path with file we just created
            input_path = INPUT_AS_VCF
    # endregion

    # region: split & annotate VCF
    if not to_path(ANNOTATED_MT).exists():
        # run annotation, default values from RefData

        siteonly_vcf_path = to_path(output_path('siteonly.vcf.gz', 'tmp'))
        input_vcf_in_batch = get_batch().read_input_group(
            **{
                'vcf.gz': str(input_path),
                'vcf.gz.tbi': str(input_path) + '.tbi',
            }
        )
        logging.info(input_path)
        sites_job, _siteonly_resource_group = add_make_sitesonly_job(
            b=get_batch(),
            input_vcf=input_vcf_in_batch,
            output_vcf_path=siteonly_vcf_path,
            storage_gb=get_config()['workflow'].get('vcf_size_in_gb', 150) + 10,
        )

        # Hack to get around previous bug in cpg-workflows that was missing units on storage reqs
        # TODO, check if this is still relevant.
        if sites_job:
            sites_job.storage(f"{get_config()['workflow'].get('vcf_size_in_gb', 150) + 10}Gi")


        # set the job dependency and cycle the 'prior' job
        if prior_job:
            sites_job.depends_on(prior_job)

        vep_ht_tmp = output_path('vep_annotations.ht', 'tmp')

        # generate the jobs which run VEP & collect the results
        vep_jobs = add_vep_jobs(
            b=get_batch(),
            input_siteonly_vcf_path=siteonly_vcf_path,
            tmp_prefix=to_path(output_path('vep_temp', 'tmp')),
            scatter_count=get_config()['workflow'].get('scatter_count', 50),
            out_path=to_path(vep_ht_tmp),
        )

        # assign sites-only job as an annotation dependency
        for job in vep_jobs:
            job.depends_on(sites_job)

        # Apply the HT of annotations to the VCF, save as MT
        anno_job = annotate_cohort_jobs(
            b=get_batch(),
            vcf_path=to_path(input_path),
            vep_ht_path=to_path(vep_ht_tmp),
            out_mt_path=to_path(ANNOTATED_MT),
            checkpoint_prefix=to_path(output_path('annotation_temp', 'tmp')),
            depends_on=vep_jobs,
            use_dataproc=False,
        )
        prior_job = anno_job[-1]

    # endregion

    #  region: query panelapp
    if not to_path(PANELAPP_JSON_OUT).exists():
        prior_job = handle_panelapp_job(
            participant_panels=participant_panels,
            prior_job=prior_job,
            previous=previous,
        )
    # endregion

    # read the ped file into the Batch
    pedigree_in_batch = get_batch().read_input(pedigree)

    # region: hail categorisation
    if not to_path(HAIL_VCF_OUT).exists():
        logging.info(f"The Labelled VCF {HAIL_VCF_OUT!r} doesn't exist; regenerating")
        prior_job = handle_hail_filtering(
            prior_job=prior_job, plink_file=pedigree_in_batch, clinvar=clinvar_table
        )
    # endregion

    # read VCF into the batch as a local file
    labelled_vcf_in_batch = (
        get_batch().read_input_group(vcf=HAIL_VCF_OUT, tbi=HAIL_VCF_OUT + '.tbi').vcf
    )

    # region: singleton decisions
    # if singleton PED supplied, also run as singletons w/separate outputs
    analysis_rounds = [(pedigree_in_batch, 'default')]
    if singletons and to_path(singletons).exists():
        to_path(singletons).copy(
            output_path(f'singletons_{EXECUTION_TIME}.fam', 'analysis')
        )
        pedigree_singletons = get_batch().read_input(singletons)
        analysis_rounds.append((pedigree_singletons, 'singletons'))
    # endregion

    # region: run results job
    # pointing this analysis at the updated config file, including input metadata
    for relationships, analysis_index in analysis_rounds:
        logging.info(f'running analysis in {analysis_index} mode')
        handle_results_job(
            labelled_vcf=labelled_vcf_in_batch,
            pedigree=relationships,
            input_path=input_path,
            output_dict=output_dict[analysis_index],
            prior_job=prior_job,
            participant_panels=participant_panels,
        )
    # endregion

    # region: copy data out
    # if we ran with per-participant panel data, copy to output folder
    # include datetime to differentiate output files and prevent clashes
    if participant_panels:
        to_path(participant_panels).copy(
            output_path(f'pid_to_panels_{EXECUTION_TIME}.json', 'analysis')
        )

    # write pedigree content to the output folder
    to_path(pedigree).copy(output_path(f'pedigree_{EXECUTION_TIME}.fam', 'analysis'))
    # endregion

    get_batch().run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )
    logging.info(
        r"""Welcome To The
          ___  _____ ______
         / _ \|_   _|| ___ \
        / /_\ \ | |  | |_/ /
        |  _  | | |  |  __/
        | | | |_| |_ | |
        \_| |_/\___/ \_|
        """
    )

    parser = ArgumentParser()
    parser.add_argument('-i', help='variant data to analyse', required=True)
    parser.add_argument('--pedigree', help='in Plink format', required=True)
    parser.add_argument('--singletons', help='singletons in Plink format')
    parser.add_argument('--participant_panels', help='per-participant panel details')
    parser.add_argument(
        '--previous',
        help='JSON file containing Gene Panel details from a prior run',
        default=None,
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
        participant_panels=args.participant_panels,
        previous=args.previous,
        skip_annotation=args.skip_annotation,
        singletons=args.singletons,
    )
