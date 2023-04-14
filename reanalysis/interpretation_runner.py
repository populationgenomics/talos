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
import sys
from argparse import ArgumentParser
from datetime import datetime
from os.path import join

from hailtop.batch.job import BashJob, Job

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import (
    authenticate_cloud_credentials_in_job,
    copy_common_env,
    dataset_path,
    output_path,
)
from cpg_workflows.batch import get_batch
from cpg_workflows.jobs.joint_genotyping import add_make_sitesonly_job
from cpg_workflows.jobs.seqr_loader import annotate_cohort_jobs
from cpg_workflows.jobs.vep import add_vep_jobs

from reanalysis import (
    hail_filter_and_label,
    html_builder,
    metamist_registration,
    mt_to_vcf,
    query_panelapp,
    summarise_clinvar_entries,
    validate_categories,
)
from reanalysis.utils import FileTypes, identify_file_type


# region: CONSTANTS
# exact time that this run occurred
EXECUTION_TIME = f'{datetime.now():%Y-%m-%d_%H:%M}'

# static paths to write outputs
ANNOTATED_MT = output_path('annotated_variants.mt')
HAIL_VCF_OUT = output_path('hail_categorised.vcf.bgz', 'analysis')
INPUT_AS_VCF = output_path('prior_to_annotation.vcf.bgz')
PANELAPP_JSON_OUT = output_path('panelapp_data.json', 'analysis')
# endregion


def set_job_resources(
    job: Job,
    prior_job: Job | None = None,
    memory: str = 'standard',
    storage: str = '20GiB',
):
    """
    apply resources to the job

    Args:
        job (Job): the job to set resources on
        prior_job (Job): the job to depend on (or None)
        memory (str): lowmem/standard/highmem
        storage (str): storage setting to use
    """
    # apply all settings to this job
    job.cpu(2).image(get_config()['workflow']['driver_image']).memory(memory).storage(
        storage
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

    if clinvar_table := get_config()['workflow'].get('forced_clinvar'):
        logging.info(f'This run will be using the Clinvar data in {clinvar_table}')
        if to_path(clinvar_table).exists():
            return None, clinvar_table

    # is it time to re-process clinvar?
    clinvar_prefix = dataset_path(
        f'clinvar_summaries/{datetime.now().strftime("%Y_%m")}'
    )
    clinvar_summary = join(clinvar_prefix, 'clinvar.ht')

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
    get_batch().write_output(bash_job.subs, join(clinvar_prefix, sub_file))
    get_batch().write_output(bash_job.vars, join(clinvar_prefix, var_file))

    # create a job to run the summary
    summarise = get_batch().new_job(name='summarise clinvar')
    set_job_resources(summarise, prior_job=bash_job)
    summarise.command(
        f'python3 {summarise_clinvar_entries.__file__} '
        f'-s {bash_job.subs} '
        f'-v {bash_job.vars} '
        f'-o {join(clinvar_prefix, "clinvar.ht")}'
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

    cmd = f'python3 {mt_to_vcf.__file__} --input {input_file} --output {INPUT_AS_VCF}'

    logging.info(f'Command used to convert MT: {cmd}')
    job.command(cmd)
    return job


def handle_panelapp_job(
    participant_panels: str | None = None, prior_job: Job | None = None
) -> Job:
    """
    creates and runs the panelapp query job

    Args:
        participant_panels (str):
        prior_job ():

    Returns:
        the Job, which other parts of the workflow may become dependent on
    """

    panelapp_job = get_batch().new_job(name='query panelapp')
    set_job_resources(panelapp_job, prior_job=prior_job)
    copy_common_env(panelapp_job)

    query_cmd = f'python3 {query_panelapp.__file__} --out_path {PANELAPP_JSON_OUT} '

    if participant_panels is not None:
        query_cmd += f'--panels {participant_panels} '

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
    set_job_resources(labelling_job, prior_job=prior_job, memory='32GiB')
    labelling_command = (
        f'python3 {hail_filter_and_label.__file__} '
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
    output: str,
    prior_job: Job | None = None,
    participant_panels: str | None = None,
):
    """
    one container to run the MOI checks, and the presentation

    Args:
        labelled_vcf (str): path to the VCF created by Hail runtime
        pedigree (str): path to the pedigree file
        input_path (str): path to the input file, logged in metadata
        output (str): path to JSON file to write
        prior_job (Job): to depend on, or None
        participant_panels (str): Optional, path to pheno-matched panels
    """

    results_job = get_batch().new_job(name='MOI tests')
    set_job_resources(results_job, prior_job=prior_job)

    gene_filter_files = (
        f'--participant_panels {participant_panels} ' if participant_panels else ''
    )

    results_command = (
        f'python3 {validate_categories.__file__} '
        f'--labelled_vcf {labelled_vcf} '
        f'--panelapp {PANELAPP_JSON_OUT} '
        f'--pedigree {pedigree} '
        f'--out_json {output} '
        f'--input_path {input_path} '
        f'{gene_filter_files}'
    )
    logging.info(f'Results command: {results_command}')
    results_job.command(results_command)
    return results_job


def handle_result_presentation_job(
    prior_job: Job | None = None, **kwargs
) -> Job | None:
    """
    run the presentation element
    allow for selection of the presentation script and its arguments

    Note: The model here is to allow other non-CPG sites/users to
          implement their own presentation scripts, and to allow for the
          same method to be shared by all users. The contract here is that
          the presentation script must one or more named arguments, and
          the argument name must match the name passed to this method,
          i.e. `--panelapp {kwargs['panelapp']} `

          The contract also requires that new presentation scripts must
          not inactivate others, i.e. a runtime configuration setting
          should allow a user to select from any of the available scripts
          based on a config parameter.

    This isn't a super slick implementation, as there are no other users
    of this method, but it's a start. Alternative scripts will have to be
    created in code, at which point the script and this little mapping will
    both have to be updated.

    kwargs currently in use:
        - results: the JSON of results created by validate_categories.py
        - panelapp: the JSON of panelapp data
        - pedigree: the pedigree file (file accessible within the batch)
        - output: the output file path

    Args:
        prior_job (Job): used in workflow dependency setting
        kwargs (): key-value arguments for presentation script
    Returns:
        The associated job
    """

    # if a new script is added, it needs to be registered here to become usable
    scripts_and_inputs = {
        'cpg': (html_builder.__file__, ['results', 'panelapp', 'pedigree', 'output'])
    }

    output_mode = get_config()['workflow'].get('presentation', 'cpg')

    # if we don't have a valid method, return the prior job
    if output_mode not in scripts_and_inputs:
        logging.warning(f'Invalid presentation mode: {output_mode}')
        return prior_job

    # extract the required inputs for the selected script
    script, required_inputs = scripts_and_inputs[output_mode]

    # set the panelapp global variable in kwargs
    kwargs['panelapp'] = PANELAPP_JSON_OUT

    display = get_batch().new_job(name='Result Presentation')
    set_job_resources(display, prior_job=prior_job)

    # assemble command from relevant variables
    script_params = ' '.join(
        [f'--{named_input} {kwargs[named_input]} ' for named_input in required_inputs]
    )
    html_command = f'python3 {script} {script_params}'

    logging.info(f'HTML generation command: {html_command}')
    display.command(html_command)
    return display


def handle_registration_jobs(
    files: list[str], registry: str, pedigree: str, prior_job: Job | None = None
):
    """
    Take a list of files and register them using the defined method.
    This registration is within a metadata DB, used to track analysis
    products, and the samples they correspond to.

    Note: Similar contract to the one as defined above in
          handle_result_presentation_job - the `registrars` mapping
          should contain all valid registration scripts, and an ID
          for each. At runtime the user should be able to select any
          registered scripts using a config parameter. If an invalid
          registrar is selected, nothing will be done.

    Args:
        files (list[str]): all files to register from this analysis
        registry (str): registration service to use
        pedigree (str): path to a pedigree file (copied into batch)
        prior_job (Job): set workflow dependency if required
    """

    # dictionary with all known registration services
    registrars = {'metamist': metamist_registration.__file__}

    # if we don't have a valid method, return the prior job
    if registry not in registrars:
        logging.error(f'Invalid registration mode: {registry}')
        return

    # create a new job that will run even if the rest of the workflow fails
    registration_job = get_batch().new_job(name='register_results')
    set_job_resources(registration_job, prior_job=prior_job)
    registration_job.always_run(True)

    metadata_command = (
        f'python3 {registrars[registry]} --pedigree {pedigree} {" ".join(files)}'
    )

    logging.info(f'Metadata registration command: {metadata_command}')
    registration_job.command(metadata_command)


def main(
    input_path: str,
    pedigree: str,
    participant_panels: str | None,
    singletons: bool = False,
    skip_annotation: bool = False,
):
    """
    main method, which runs the full reanalysis process

    Args:
        input_path (str): path to the VCF/MT
        pedigree (str): family file for this analysis
        participant_panels (str): file containing panels-per-family (optional)
        singletons (bool): run as Singletons (with appropriate output paths)
        skip_annotation (bool): if the input is annotated, don't re-run
    """

    assert to_path(
        input_path
    ).exists(), f'The provided path {input_path!r} does not exist or is inaccessible'

    logging.info('Starting the reanalysis batch')

    # region: output files lookup
    # separate paths for familial and singleton analysis
    if singletons:
        assert (
            'singleton' in get_config()['workflow']['output_prefix']
        ), 'To keep singletons separate, include "singleton" in the file path'

    # modify output paths depending on analysis type
    output_dict = {
        'web_html': output_path(
            f'{"singleton" if singletons else "summary"}_output.html', 'web'
        ),
        'results': output_path(
            f'{"singleton" if singletons else "summary"}_results.json', 'analysis'
        ),
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
        output_dict['annotated_mt'] = ANNOTATED_MT
        prior_job = anno_job[-1]
    # endregion

    #  region: query panelapp
    if not to_path(PANELAPP_JSON_OUT).exists():
        prior_job = handle_panelapp_job(
            participant_panels=participant_panels, prior_job=prior_job
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
        output_dict['hail_vcf'] = HAIL_VCF_OUT
    # endregion

    # read VCF into the batch as a local file
    labelled_vcf_in_batch = (
        get_batch().read_input_group(vcf=HAIL_VCF_OUT, tbi=HAIL_VCF_OUT + '.tbi').vcf
    )

    # region: run results job
    # pointing this analysis at the updated config file, including input metadata
    prior_job = handle_results_job(
        labelled_vcf=labelled_vcf_in_batch,
        pedigree=pedigree_in_batch,
        input_path=input_path,
        output=output_dict['results'],
        prior_job=prior_job,
        participant_panels=participant_panels,
    )
    prior_job = handle_result_presentation_job(
        prior_job=prior_job,
        pedigree=pedigree_in_batch,
        output=output_dict['web_html'],
        results=output_dict['results'],
    )
    prior_job.always_run()
    # endregion

    # region: output registration job
    # register the output files in metamist if required
    if registry := get_config()['workflow'].get('register'):
        logging.info(f'Metadata registration will be done using {registry}')
        handle_registration_jobs(
            files=sorted(output_dict.values()),
            pedigree=pedigree_in_batch,
            registry=registry,
            prior_job=prior_job,
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
    parser.add_argument('--participant_panels', help='per-participant panel details')
    parser.add_argument(
        '--singletons',
        help='boolean, set if this run is a singleton pedigree',
        action='store_true',
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
        skip_annotation=args.skip_annotation,
        singletons=args.singletons,
    )
