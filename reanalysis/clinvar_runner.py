#!/usr/bin/env python3


"""
Entrypoint for clinvar summary generation
"""

from datetime import datetime

import click
from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import (
    authenticate_cloud_credentials_in_job,
    output_path,
    query_command,
)
from cpg_workflows.batch import get_batch

from reanalysis import clinvar_by_codon, summarise_clinvar_entries, seqr_loader
from reanalysis.vep_jobs import add_vep_jobs


@click.command
@click.option('--ht_out', help='Path to write the Hail table to')
@click.option('--date', help='date cut-off, optional', default=None)
def main(ht_out: str, date: str | None = None):
    """
    run the clinvar summary, output to defined path
    Args:
        ht_out ():
        date ():
    """

    clinvar_table_path = to_path(ht_out)
    clinvar_folder = clinvar_table_path.parent

    # create a bash job to copy data from remote
    bash_job = get_batch().new_bash_job(name='copy clinvar files to local')
    bash_job.image(get_config()['workflow']['driver_image'])

    directory = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/'
    sub_file = 'submission_summary.txt.gz'
    var_file = 'variant_summary.txt.gz'

    bash_job.command(
        (
            f'wget -q {directory}{sub_file} -O {bash_job.subs} && '
            f'wget -q {directory}{var_file} -O {bash_job.vars}'
        )
    )

    # construct a date String to identify the origina date for these files
    today = datetime.now().strftime('%y-%m-%d')

    # write output files date-specific
    get_batch().write_output(bash_job.subs, str(clinvar_folder / f'{today}_{sub_file}'))
    get_batch().write_output(bash_job.vars, str(clinvar_folder / f'{today}_{var_file}'))

    # create a space for the SNV VCF
    snv_vcf = clinvar_folder / 'pathogenic_snv.vcf.bgz'

    # region: run the summarise_clinvar_entries script
    summarise = get_batch().new_job(name='summarise clinvar')
    prior_job = summarise

    summarise.cpu(2).image(get_config()['workflow']['driver_image']).storage('20G')
    authenticate_cloud_credentials_in_job(summarise)
    command_options = (
        f'-s {bash_job.subs} '
        f'-v {bash_job.vars} '
        f'-o {ht_out} '
        f'--path_snv {snv_vcf} '
    )
    if date:
        command_options += f' -d {date}'
    summarise.command(f'python3 {summarise_clinvar_entries.__file__} {command_options}')
    # endregion

    # region: annotate the SNV VCF with VEP
    vep_ht_tmp = to_path(output_path('vep_annotations.ht', 'tmp'))
    annotated_clinvar = to_path(clinvar_folder / 'annotated_clinvar.mt')

    # generate the jobs which run VEP & collect the results
    vep_jobs = add_vep_jobs(
        b=get_batch(),
        input_siteonly_vcf_path=snv_vcf,
        tmp_prefix=to_path(output_path('vep_temp', 'tmp')),
        scatter_count=5,
        out_path=vep_ht_tmp,
    )

    # add Clinvar job as an annotation dependency
    if prior_job and vep_jobs:
        for job in vep_jobs:
            job.depends_on(prior_job)
        prior_job = vep_jobs[-1]

    j = get_batch().new_job(f'annotate cohort')
    j.image(get_config()['workflow']['driver_image'])
    j.command(
        query_command(
            seqr_loader,
            seqr_loader.annotate_cohort.__name__,
            str(snv_vcf),
            str(annotated_clinvar),
            str(vep_ht_tmp),
            output_path('annotation_temp', 'tmp'),
            setup_gcp=True,
        )
    )
    if prior_job:
        j.depends_on(prior_job)
    prior_job = j
    # endregion

    # region: run the clinvar_by_codon script

    pm5_table = to_path(clinvar_folder / 'pm5_table.mt')
    clinvar_by_codon_job = get_batch().new_job(name='clinvar_by_codon')
    clinvar_by_codon_job.image(get_config()['workflow']['driver_image']).cpu(2).storage(
        '20G'
    )
    authenticate_cloud_credentials_in_job(clinvar_by_codon_job)
    clinvar_by_codon_job.command(
        f'python3 {clinvar_by_codon.__file__} '
        f'--mt {annotated_clinvar} '
        f'--write_path {pm5_table}'
    )
    clinvar_by_codon_job.depends_on(prior_job)
    # endregion
    get_batch().run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
