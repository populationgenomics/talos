#!/usr/bin/env python3


"""
Entrypoint for clinvar summary generation
"""

from datetime import datetime
from os.path import join

from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import authenticate_cloud_credentials_in_job, get_batch, query_command

from reanalysis import clinvar_by_codon, seqr_loader, summarise_clinvar_entries
from reanalysis.static_values import get_logger
from reanalysis.vep_jobs import add_vep_jobs


def generate_clinvar_table(cloud_folder: Path, clinvar_outputs: str):
    """
    set up the job that does de novo clinvar summary

    Args:
        cloud_folder (Path): folder for this analysis
        clinvar_outputs (str): prefix for writing new files/dirs
    """

    bash_job = get_batch().new_bash_job(name='copy clinvar files to local')
    bash_job.image(config_retrieve(['workflow', 'driver_image']))

    directory = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/'
    sub_file = 'submission_summary.txt.gz'
    var_file = 'variant_summary.txt.gz'

    bash_job.command(
        f'wget -q {directory}{sub_file} -O {bash_job.subs} && wget -q {directory}{var_file} -O {bash_job.vars}',
    )

    # write output files date-specific
    get_batch().write_output(bash_job.subs, str(cloud_folder / sub_file))
    get_batch().write_output(bash_job.vars, str(cloud_folder / var_file))

    # region: run the summarise_clinvar_entries script
    summarise = get_batch().new_job(name='summarise clinvar')
    summarise.depends_on(bash_job)

    summarise.cpu(2).image(config_retrieve(['workflow', 'driver_image'])).storage('20G')
    authenticate_cloud_credentials_in_job(summarise)
    summarise.command(
        f'python3 {summarise_clinvar_entries.__file__} '
        f'-s {bash_job.subs} '
        f'-v {bash_job.vars} '
        f'-o {clinvar_outputs}',
    )

    return summarise


def generate_annotated_data(annotation_out: Path, snv_vcf: str, tmp_path: Path, dependency: Job | None = None) -> Job:
    """
    if the annotated data Table doesn't exist, generate it

    Args:
        annotation_out (str): MT path to create
        snv_vcf (str): path to a VCF file
        tmp_path (Path): path to a temporary folder
        dependency (Job | None): optional job dependency

    Returns:
        The Job for future dependency setting
    """

    vep_ht_tmp = tmp_path / 'vep_annotations.ht'

    # generate the jobs which run VEP & collect the results
    vep_jobs = add_vep_jobs(
        b=get_batch(),
        input_siteonly_vcf_path=snv_vcf,
        tmp_prefix=tmp_path / 'vep_temp',
        scatter_count=50,
        out_path=vep_ht_tmp,
    )

    # add Clinvar job as an annotation dependency
    # update dependency job if necessary
    if vep_jobs:
        if dependency:
            for job in vep_jobs:
                job.depends_on(dependency)
        dependency = vep_jobs[-1]

    j = get_batch().new_job('annotate cohort')
    j.image(config_retrieve(['workflow', 'driver_image']))

    # run seqr_loader, only applying VEP annotations
    j.command(
        query_command(
            seqr_loader,
            seqr_loader.annotate_cohort.__name__,
            str(snv_vcf),
            str(annotation_out),
            str(vep_ht_tmp),
            str(tmp_path / 'annotation_temp'),
            True,
            setup_gcp=True,
        ),
    )
    if dependency:
        j.depends_on(dependency)
    return j


def main():
    """
    run the clinvar summary, output to common path
    """

    cloud_folder = to_path(
        join(config_retrieve(['storage', 'common', 'analysis']), 'aip_clinvar', datetime.now().strftime('%y-%m')),
    )

    # clinvar VCF, decisions, annotated VCF, and PM5
    clinvar_output_path = join(str(cloud_folder), 'clinvar_decisions')
    clinvar_ht = f'{clinvar_output_path}.ht'
    snv_vcf = f'{clinvar_output_path}.vcf.bgz'
    clinvar_pm5_path = cloud_folder / 'clinvar_pm5.ht'
    annotated_clinvar = cloud_folder / 'annotated_clinvar.mt'

    # check if we can just quit already
    if all(this_path.exists() for this_path in [annotated_clinvar, clinvar_ht, clinvar_pm5_path]):
        get_logger().info('Clinvar data already exists, exiting')
        return

    temp_path = to_path(
        join(config_retrieve(['storage', 'common', 'tmp']), 'aip_clinvar', datetime.now().strftime('%y-%m')),
    )

    dependency = None

    # generate a new round of clinvar decisions
    if not all(to_path(output).exists() for output in [clinvar_ht, snv_vcf]):
        dependency = generate_clinvar_table(cloud_folder, clinvar_output_path)

    # create the annotation job(s)
    if not annotated_clinvar.exists():
        dependency = generate_annotated_data(annotated_clinvar, snv_vcf, temp_path, dependency)

    # region: run the clinvar_by_codon script
    if not clinvar_pm5_path.exists():
        clinvar_by_codon_job = get_batch().new_job(name='clinvar_by_codon')
        clinvar_by_codon_job.image(config_retrieve(['workflow', 'driver_image'])).cpu(2).storage('20G')
        authenticate_cloud_credentials_in_job(clinvar_by_codon_job)
        clinvar_by_codon_job.command(
            f'python3 {clinvar_by_codon.__file__} --mt_path {annotated_clinvar} --write_path {clinvar_pm5_path}',
        )
        if dependency:
            clinvar_by_codon_job.depends_on(dependency)
    # endregion

    get_batch().run(wait=False)


if __name__ == '__main__':
    main()
